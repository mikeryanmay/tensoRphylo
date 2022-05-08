/*
 * EigenKernels.cpp
 *
 *  Created on: Sep 2, 2019
 *      Author: xaviermeyer
 */

#include "EigenKernels.h"

#include "../EigenUtils.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/StateTypes/Utils.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/Scheduler/DAG/NodeDAG.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/Approximator/Specialized/BaseDense.h"
#include "Likelihood/Approximator/Specialized/AdaptivePairedUS.h"
#include "Likelihood/Approximator/Specialized/DenseUnobserved.h"
#include "Likelihood/Approximator/Specialized/SpecializedApproxManager.h"

#include <Eigen/Core>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Branchwise {

EigenKernels::EigenKernels(Conditions::conditionalProbability_t aConditionType,
							Phylogeny::Data::ContainerSharedPtr aPtrData,
							 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager,
							SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
							Tensor::ContainerSharedPtr aPtrTensorsContainer,
							bool &aCondCompatibility) :
									condCompatibility(aCondCompatibility),
									conditionType(aConditionType),
									ptrData(aPtrData),
									ptrSApproxManager(aPtrSApproxManager),
									ptrSynchEventContainer(aPtrSynchEventContainer),
									ptrTensorsContainer(aPtrTensorsContainer) {
}

EigenKernels::~EigenKernels() {
}


void EigenKernels::setInitialCondition(Likelihood::Scheduler::DAG::NodeDAG *node) {

	// Get initial mass sampling event
	const Eigen::VectorXd& vecRhoOf0 = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(0.);

	// Get the data and allocate memory in the state type
	assert(!node->getEvent()->getNodes().empty() && node->getEvent()->getNodes().front()->isExtant());
	size_t idTaxa = node->getEvent()->getNodes().front()->getTaxaId();
	const Eigen::VectorXd &probs = ptrData->getProbForTaxaIdThreadSafe(idTaxa);

	Likelihood::StateType::Vector::EigenState *pState = new Likelihood::StateType::Vector::EigenState();
	pState->allocateVecProb();
	Eigen::VectorXd& p = pState->getStateProb();
	p = probs.cwiseProduct(vecRhoOf0);

	// Attach the state to the dag node
	pState->setEdgeMapping(node->getEdgeIdMapping());
	assert(node->getProbStates().size() == 0);
	node->attachState(pState);

}

void EigenKernels::setInitalExtinctNodeCondition(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {

	// Get the observed data
	assert(!node->getEvent()->getNodes().empty() && node->getEvent()->getNodes().front()->isExtinct());
	size_t idTaxa = node->getEvent()->getNodes().front()->getTaxaId();
	const Eigen::VectorXd &probs = ptrData->getProbForTaxaIdThreadSafe(idTaxa);

	// Allocate a new vector of prob for the new fossil sample/parent edge
	Likelihood::StateType::Vector::EigenState *pState = new Likelihood::StateType::Vector::EigenState();
	pState->allocateVecProb();

	// Get parameters
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	// Get the prob vector and initialize it
	Eigen::VectorXd& p = pState->getStateProb();
	Eigen::VectorXd u;
	ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(t,u);
	p = probs.cwiseProduct(phi.cwiseProduct(u) + delta);

	// Attach the state to the dag node
	pState->setEdgeMapping(node->getEdgeIdMapping());
	node->attachState(pState);

}

void EigenKernels::computeAsynchSampling(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {
	assert(!node->getEvent()->getNodes().empty() && node->getEvent()->getNodes().front()->isSampledAncestor());

	// Get parameters
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);

	// Get the observed data
	size_t idTaxa = node->getEvent()->getNodes().front()->getTaxaId();
	const Eigen::VectorXd &probs = ptrData->getProbForTaxaIdThreadSafe(idTaxa);

	// Get the vector as a function of the before edge ID
	assert(node->getProbStates().size() == 1);
	Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
	// If the child node of this node IS NOT a ghost node: we compute as usual the lik
	// If the child node of this node IS a ghost node: we force the state probabilities to u
	Eigen::VectorXd tmpProb;
	if(!node->getEvent()->getNodes().front()->getChildrenNodes()[0]->isGhostNode()) { // If its node a ghost node, we just compute as usual
		tmpProb = p;
	} else {
		ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(t,tmpProb);
	}
	p = phi.cwiseProduct(tmpProb.cwiseProduct(probs));
	rescaleProbabilities(node);
}

void EigenKernels::computeAsynchSpeciation(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {
	assert(!node->getEvent()->getNodes().empty() && node->getEvent()->getNodes().front()->isSpeciationNode());
	assert(node->getProbStates().size() == 2);

	// Get children vec prob
	Eigen::VectorXd& obsStatesProbL = node->getProbStates().front()->getStateProb();
	Eigen::VectorXd& obsStatesProbR = node->getProbStates().back()->getStateProb();

	// Get parameters and tensors
	Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);

	// Do the contractions
	if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
		Eigen::VectorXd resContraction = obsStatesProbL.cwiseProduct(obsStatesProbR);
		obsStatesProbL = resContraction.cwiseProduct(lambda);
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		Eigen::VectorXd resContraction = computeTensorContractionVector(omega, obsStatesProbL, obsStatesProbR);
		obsStatesProbL = resContraction.cwiseProduct(lambda);
	} else {
		Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
		Eigen::VectorXd resContraction = computeTensorContractionVector(omega, obsStatesProbL, obsStatesProbR);
		obsStatesProbL = resContraction.cwiseProduct(lambda);
	}

	// Get scaling factors:
	double scalingProbL = node->getProbStates().front()->getScaling();
	double scalingProbR = node->getProbStates().back()->getScaling();
	double newScaling = scalingProbL+scalingProbR;
	newScaling += Likelihood::StateType::rescaleProbabilityVector(obsStatesProbL);
	// Set new scaling factor
	node->getProbStates().front()->setScaling(newScaling);

	// remove and delete the state of the right branch (no longer needed)
	node->getProbStates().back()->removeVecProb(); // This is already done in the dtor - to be safe
	delete node->getProbStates().back();
	node->getProbStates().pop_back();
}

void EigenKernels::computeMassSpeciation(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {

	// Get mass speciation event prob
	const Eigen::VectorXd& vecUpsilonOfT = ptrSynchEventContainer->getPtrMassSpeciation()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusUpsilonOfT = Eigen::VectorXd::Ones(vecUpsilonOfT.size()) - vecUpsilonOfT;

	if(node->getProbStates().size() == 1) { // its an edge
		Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
		Eigen::VectorXd u;
		ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(t,u);

		// Do the contractions
		if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
			Eigen::VectorXd resContraction = p.cwiseProduct(u);
			p = p.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContraction);
		} else if(ptrTensorsContainer->getOmega()->isSparse()) {
			Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
			// Do the contractions
			Eigen::VectorXd resContraction = computeTensorContractionVector(omega, p, u);
			p = p.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContraction);
		} else {
			Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
			// Do the contractions
			Eigen::VectorXd resContraction = computeTensorContractionVector(omega, p, u);
			p = p.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContraction);
		}
		rescaleProbabilities(node);
	} else if(node->getProbStates().size() == 2) { // its an speciation event
		Eigen::VectorXd& pL = node->getProbStates().front()->getStateProb();
		Eigen::VectorXd& pR = node->getProbStates().back()->getStateProb();

		if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
			Eigen::VectorXd resContraction = pL.cwiseProduct(pR);
			pL = vecUpsilonOfT.cwiseProduct(resContraction);
		} else if(ptrTensorsContainer->getOmega()->isSparse()) {
			Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
			// Do the contraction
			Eigen::VectorXd resContraction = computeTensorContractionVector(omega, pL, pR);
			// Storing in pL, pR will have to be removed afterward
			pL = vecUpsilonOfT.cwiseProduct(resContraction);
		} else {
			Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
			// Do the contraction
			Eigen::VectorXd resContraction = computeTensorContractionVector(omega, pL, pR);
			// Storing in pL, pR will have to be removed afterward
			pL = vecUpsilonOfT.cwiseProduct(resContraction);
		}

		// Get scaling factors:
		double scalingProbL = node->getProbStates().front()->getScaling();
		double scalingProbR = node->getProbStates().back()->getScaling();
		double newScaling = scalingProbL+scalingProbR;
		newScaling += Likelihood::StateType::rescaleProbabilityVector(pL);
		// Set new scaling factor
		node->getProbStates().front()->setScaling(newScaling);

		// remove and delete the state of the right branch (no longer needed)
		node->getProbStates().back()->removeVecProb(); // This is already done in the dtor - to be safe
		delete node->getProbStates().back();
		node->getProbStates().pop_back();

	} else { // its wrong
		assert(false && "Mass speciation should be affecting edges or speciation nodes.");
	}
}

void EigenKernels::computeMassExtinction(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {

	assert(node->getEvent()->getNodes().empty() && node->getProbStates().size() == 1);

	if(!ptrSynchEventContainer->getPtrMassExtinction()->areStateChangeInvolved()) { // No state change prob
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// We compute for edges -  some work to do here for openmp implementation
		Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
		p = vecOneMinusGammaOfT.cwiseProduct(p);
	} else {
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		const Eigen::MatrixXd& matZetaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getStateChangeProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// We compute for edges -  some work to do here for openmp implementation
		Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
		p = vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*p);
	}
	rescaleProbabilities(node);

}

void EigenKernels::computeMassSamplingEvent(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {

	// Sampling probability
	const Eigen::VectorXd& vecRhoOfT = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusRhoOfT = Eigen::VectorXd::Ones(vecRhoOfT.size()) - vecRhoOfT;

	if(node->getEvent()->getNodes().empty()) { // Lineage
		// Observed/unobserved probabilities
		assert(node->getProbStates().size() == 1);
		Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();

		p = vecOneMinusRhoOfT.cwiseProduct(p);
	} else if (node->getEvent()->getNodes().size() == 1 && node->getEvent()->getNodes().front()->isExtinct()) {
		assert(node->getProbStates().size() == 0);

		// Get the observed data
		size_t idTaxa = node->getEvent()->getNodes().front()->getTaxaId();
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaIdThreadSafe(idTaxa);

		// Get unobserved prob
		Eigen::VectorXd u;
		ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(t,u);

		// Allocate a new vector of prob for the new fossil sample/parent edge
		Likelihood::StateType::Vector::EigenState *pState = new Likelihood::StateType::Vector::EigenState();
		pState->allocateVecProb();

		// Get the prob vector and initialize it
		Eigen::VectorXd& p = pState->getStateProb();
		p = probs.cwiseProduct(vecRhoOfT.cwiseProduct(u));

		pState->setEdgeMapping(node->getEdgeIdMapping());
		node->attachState(pState);

	} else if (node->getEvent()->getNodes().size() == 1 && node->getEvent()->getNodes().front()->isSampledAncestor()) {
		assert(node->getProbStates().size() == 1);
		// Get the observed data
		size_t idTaxa = node->getEvent()->getNodes().front()->getTaxaId();
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaIdThreadSafe(idTaxa);

		// Get the vector as a function of the before edge ID
		Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
		Eigen::VectorXd tmpP = p; // Copy
		p = probs.cwiseProduct(vecRhoOfT.cwiseProduct(tmpP));

	} else {
		assert(false && "This event does not match the possible scenarii for mass sampling events.");
	}
	rescaleProbabilities(node);
}

void EigenKernels::computeMassDestrSamplingEvent(double t, Likelihood::Scheduler::DAG::NodeDAG *node) {

	// Sampling probability
	const Eigen::VectorXd& vecXiOfT = ptrSynchEventContainer->getPtrMassDestrSampling()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusXiOfT = Eigen::VectorXd::Ones(vecXiOfT.size()) - vecXiOfT;

	if(node->getEvent()->getNodes().empty()) {
		assert(node->getProbStates().size() == 1);
		Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
		p = vecOneMinusXiOfT.cwiseProduct(p);
	} else if (node->getEvent()->getNodes().size() == 1 && node->getEvent()->getNodes().front()->isExtinct()) {
		assert(node->getProbStates().size() == 0);

		// Get the observed data
		size_t idTaxa = node->getEvent()->getNodes().front()->getTaxaId();
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaIdThreadSafe(idTaxa);

		// Get unobserved prob
		Eigen::VectorXd u;
		ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(t,u);

		// Allocate a new vector of prob for the new fossil sample/parent edge
		Likelihood::StateType::Vector::EigenState *pState = new Likelihood::StateType::Vector::EigenState();
		pState->allocateVecProb();

		// Get the prob vector and initialize it
		Eigen::VectorXd& p = pState->getStateProb();
		p = vecXiOfT.cwiseProduct(probs);


		pState->setEdgeMapping(node->getEdgeIdMapping());
		node->attachState(pState);
		rescaleProbabilities(node);

	} else {
		assert(false && "This event does not match the possible scenarii for mass sampling events.");
	}
}

void EigenKernels::rescaleProbabilities(Likelihood::Scheduler::DAG::NodeDAG *node) {

	assert(node->getProbStates().size() == 1);

	// Get children vec prob
	Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
	double scalingFactor = node->getProbStates().front()->getScaling();
	scalingFactor += Likelihood::StateType::rescaleProbabilityVector(p);

	// Set new scaling factor
	node->getProbStates().front()->setScaling(scalingFactor);

}

double EigenKernels::computeLogLikelihood(double t, Eigen::VectorXd &prior, Likelihood::Scheduler::DAG::NodeDAG *node) {
	assert(!node->getEvent()->getNodes().empty() && node->getEvent()->getNodes().front()->isOriginNode());
	assert(node->getProbStates().size() == 1);

	using namespace Likelihood::Approximator::Specialized;

	Eigen::VectorXd& p = node->getProbStates().front()->getStateProb();
	const Eigen::VectorXd &u = ptrSApproxManager->getSpecializedUnobserved(FULL_PROCESS)->getFinalProbability();
	assert(p.size() == prior.size());

	double unscaledPData = prior.dot(p);
	double logPData = log(unscaledPData) + node->getProbStates().front()->getScaling();
	double pSurvival = 1.0;

	if(conditionType == Likelihood::Conditions::TIME ) {
		// Do nothing
	} else if (conditionType == Likelihood::Conditions::ROOT_SURVIVAL ) {
		const Eigen::VectorXd &uHat = ptrSApproxManager->getSpecializedUnobserved(EXTANT_PROCESS)->getFinalProbability();
		Eigen::VectorXd vecOneMinusUHat = Eigen::VectorXd::Ones(uHat.size()) - uHat;
		if(condCompatibility || ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) {
			Eigen::VectorXd vecOneMinusUHat2 = vecOneMinusUHat.cwiseProduct(vecOneMinusUHat);
			pSurvival = prior.dot(vecOneMinusUHat2);
		} else {
			// Do the contractions
			if(ptrTensorsContainer->getOmega()->isSparse()) { // Sparse clado
				Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusUHat, vecOneMinusUHat);
				pSurvival = prior.dot(resContraction);
			} else {
				Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusUHat, vecOneMinusUHat);
				pSurvival = prior.dot(resContraction);
			}
		}
	} else if (conditionType == Likelihood::Conditions::ROOT_MRCA ) {
		const Eigen::VectorXd &uHat = ptrSApproxManager->getSpecializedUnobserved(EXTANT_PROCESS)->getFinalProbability();
		Eigen::VectorXd vecOneMinusUHat = Eigen::VectorXd::Ones(uHat.size()) - uHat;
		Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
		if(condCompatibility || ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) {
			Eigen::VectorXd vecOneMinusUHat2 = vecOneMinusUHat.cwiseProduct(vecOneMinusUHat);
			Eigen::VectorXd specTimesOneMinusUHat2 = vecOneMinusUHat2.cwiseProduct(lambda);
			pSurvival = prior.dot(specTimesOneMinusUHat2);
		} else {
			// Do the contractions
			if(ptrTensorsContainer->getOmega()->isSparse()) { // Sparse clado
				Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusUHat, vecOneMinusUHat);
				pSurvival = prior.dot(resContraction.cwiseProduct(lambda));
			} else {
				Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusUHat, vecOneMinusUHat);
				pSurvival = prior.dot(resContraction.cwiseProduct(lambda));
			}
		}
	} else if (conditionType == Likelihood::Conditions::ROOT_SAMPLING ) {
		Eigen::VectorXd vecOneMinusU = Eigen::VectorXd::Ones(u.size()) - u;
		if(condCompatibility || ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) {
			Eigen::VectorXd vecOneMinusU2 = vecOneMinusU.cwiseProduct(vecOneMinusU);
			pSurvival = prior.dot(vecOneMinusU2);
		} else {
			// Do the contractions
			if(ptrTensorsContainer->getOmega()->isSparse()) { // Sparse clado
				Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusU, vecOneMinusU);
				pSurvival = prior.dot(resContraction);
			} else {
				Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusU, vecOneMinusU);
				pSurvival = prior.dot(resContraction);
			}
		}
	} else if (conditionType == Likelihood::Conditions::STEM_SURVIVAL ) {
		const Eigen::VectorXd &uHat = ptrSApproxManager->getSpecializedUnobserved(EXTANT_PROCESS)->getFinalProbability();
		Eigen::VectorXd vecOneMinusUHat = Eigen::VectorXd::Ones(uHat.size()) - uHat;
		pSurvival = prior.dot(vecOneMinusUHat);
	} else if (conditionType == Likelihood::Conditions::STEM_ONE_SAMPLE ) {
		Eigen::VectorXd vecOneMinusU = Eigen::VectorXd::Ones(u.size()) - u;
		pSurvival = prior.dot(vecOneMinusU);
	} else if (conditionType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES ) {
		const Eigen::VectorXd &uHat = ptrSApproxManager->getSpecializedPairedUS(EXTANT_PROCESS)->getU();
		const Eigen::VectorXd &sHat = ptrSApproxManager->getSpecializedPairedUS(EXTANT_PROCESS)->getS();
		Eigen::VectorXd vecOneMinusUHatSHat = Eigen::VectorXd::Ones(uHat.size()) - uHat ;
		vecOneMinusUHatSHat -= sHat;
		pSurvival = prior.dot(vecOneMinusUHatSHat);
	} else if (conditionType == Likelihood::Conditions::STEM_TWO_SAMPLES ) {
		const Eigen::VectorXd &s = ptrSApproxManager->getSpecializedPairedUS(FULL_PROCESS)->getS();
		Eigen::VectorXd vecOneMinusUS = Eigen::VectorXd::Ones(u.size()) - u;
		vecOneMinusUS -= s;
		pSurvival = prior.dot(vecOneMinusUS);
	} else if (conditionType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
		Eigen::VectorXd vecOneMinusU = Eigen::VectorXd::Ones(u.size()) - u;
		Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
		if(condCompatibility || ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) {
			Eigen::VectorXd vecOneMinusU2 = vecOneMinusU.cwiseProduct(vecOneMinusU);
			pSurvival = prior.dot(vecOneMinusU2.cwiseProduct(lambda));
		} else {
			// Do the contractions
			if(ptrTensorsContainer->getOmega()->isSparse()) { // Sparse clado
				Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusU, vecOneMinusU);
				pSurvival = prior.dot(resContraction.cwiseProduct(lambda));
			} else {
				Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
				Eigen::VectorXd resContraction = computeTensorContractionVector(omega, vecOneMinusU, vecOneMinusU);
				pSurvival = prior.dot(resContraction.cwiseProduct(lambda));
			}
		}
	} else {
		assert(false && "conditionType unknown -- unsupported survival conditioning. ");
	}

	//std::cout << "pSurvival = " << pSurvival << std::endl;
	return logPData - log(pSurvival);
}

} /* namespace Branchwise */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
