/*
 * EigenKernels.cpp
 *
 *  Created on: Nov 19, 2019
 *      Author: xaviermeyer
 */

#include "EigenKernels.h"


#ifndef LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENKERNELS_DEF_H_
#define LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENKERNELS_DEF_H_

#if defined(_OPENMP)

#include "EigenUtils.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/StateTypes/Utils.h"
#include "Likelihood/Scheduler/IncScheduler.h"

#include <Eigen/Core>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace OpenMP {

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
EigenKernels<conditionalProbType, withCladoEvents>::EigenKernels(Phylogeny::Data::ContainerSharedPtr aPtrData,
							SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
							Tensor::ContainerSharedPtr aPtrTensorsContainer,
							bool &aCondCompatibility) :
									condCompatibility(aCondCompatibility),
									ptrData(aPtrData),
									ptrSynchEventContainer(aPtrSynchEventContainer),
									ptrTensorsContainer(aPtrTensorsContainer) {
}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
EigenKernels<conditionalProbType, withCladoEvents>::~EigenKernels() {
}


template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::setInitialCondition(const std::vector<PS::Node*> &extantNodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {

	size_t nNodes = extantNodes.size();

	// Get initial mass sampling event
	const Eigen::VectorXd& vecRhoOf0 = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(0.);

	// For each extant nodes
	for(size_t iN = 0; iN < nNodes; ++iN) {
		assert(extantNodes[iN]->isExtant());

		// Get the data and allocate memory in the state type
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaId(extantNodes[iN]->getTaxaId());
		size_t iVec = x.allocateVecProbForEdge(extantNodes[iN]->getEdgeToParent()->getId());
		assert(iVec == iN);

		// Get the appropriate probability vec from the state type and init it as a function of the observed state
		Eigen::Ref< Eigen::VectorXd > p = x.getObservedStateProb(iVec);
		p = probs.cwiseProduct(vecRhoOf0);
	}

	// Init the unobserved probs
	Eigen::Ref< Eigen::VectorXd > unobservedStatesProb = x.getUnobservedStateProb();
	unobservedStatesProb.setOnes();
	unobservedStatesProb -= vecRhoOf0;

	// Init what's needed
	// singleton:
	if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
		s = vecRhoOf0;
	}

	// unobserved no sampling:
	if(conditionalProbType == Conditions::ROOT_SURVIVAL || conditionalProbType == Conditions::ROOT_MRCA  ||
	   conditionalProbType == Conditions::STEM_SURVIVAL || conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();
		uHat.setOnes();
		uHat -= vecRhoOf0;
	}

	// singleton no sampling:
	if(conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > sHat = x.getSingletonNoSamplingStateProb();
		sHat = vecRhoOf0;
	}

	assert(x.size() == nNodes);
}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::setInitalExtinctNodeCondition(double t, PS::Node* extinctNode, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {
	assert(extinctNode->isExtinct());

	// Get the observed data
	const Eigen::VectorXd &probs = ptrData->getProbForTaxaId(extinctNode->getTaxaId());

	// Allocate a new vector of prob for the new fossil sample/parent edge
	size_t iVec = x.allocateVecProbForEdge(extinctNode->getEdgeToParent()->getId());

	// Get parameters
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	// Get the prob vector and initialize it
	Eigen::Ref< Eigen::VectorXd > unobservedStatesProb = x.getUnobservedStateProb();
	Eigen::Ref< Eigen::VectorXd > p = x.getObservedStateProb(iVec);
	p = probs.cwiseProduct(phi.cwiseProduct(unobservedStatesProb) + delta);

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::computeAsynchSampling(double t, PS::Node* samplingNode, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {
	assert(samplingNode->isSampledAncestor());

	// keep track of edges id
	size_t iBeforeEdge = samplingNode->getEdgesToChildren().front()->getId();
	size_t iAfterEdge = samplingNode->getEdgeToParent()->getId();

	// Get parameters
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);

	// Get the observed data
	const Eigen::VectorXd &probs = ptrData->getProbForTaxaId(samplingNode->getTaxaId());

	// Get the vector as a function of the before edge ID
	size_t iBeforeProbVector = x.findVecProbIdForEdgeId(iBeforeEdge);
	Eigen::Ref< Eigen::VectorXd > p = x.getObservedStateProb(iBeforeProbVector);
	// If the child node of this node IS NOT a ghost node: we compute as usual the lik
	// If the child node of this node IS a ghost node: we force the state probabilities to one
	Eigen::VectorXd tmpProb;
	if(!samplingNode->getChildrenNodes()[0]->isGhostNode()) { // If its node a ghost node, we just compute as usual
		tmpProb = p;
	} else {
		tmpProb = x.getUnobservedStateProb();
	}

	p = phi.cwiseProduct(tmpProb.cwiseProduct(probs));

	// Update the Edge to prob vector mapping
	x.updateVecProbForEdgeMapping(iBeforeEdge, iAfterEdge);

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::computeAsynchSpeciation(double t, PS::Node* parentNode, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {
	assert(!parentNode->isSampledAncestor() && !parentNode->isExtinct() && !parentNode->isExtant());
	assert(parentNode->getEdgesToChildren().size() == 2);

	// Identify children edges
	size_t iChildEdgeL = parentNode->getEdgesToChildren()[0]->getId();
	size_t iChildEdgeR = parentNode->getEdgesToChildren()[1]->getId();

	// Identify parent edge
	size_t iParentEdge;
	if(parentNode->getEdgeToParent() == NULL) { // Special case for root node without stem
		iParentEdge = 0;
	} else {
		iParentEdge = parentNode->getEdgeToParent()->getId();
	}

	// Get children vec prob
	Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
	size_t iChildLProbVector = x.findVecProbIdForEdgeId(iChildEdgeL);
	Eigen::Ref< Eigen::MatrixXd >::ColXpr pL = obsStatesProb.col(iChildLProbVector);
	size_t iChildRProbVector = x.findVecProbIdForEdgeId(iChildEdgeR);
	Eigen::Ref< Eigen::MatrixXd >::ColXpr pR = obsStatesProb.col(iChildRProbVector);

	// Get parameters and tensors
	Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);

	// Do the contractions
	Eigen::VectorXd resContraction = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::MatrixXd >::ColXpr, Eigen::Ref< Eigen::MatrixXd >::ColXpr >(ptrTensorsContainer->getOmega(), pL, pR, t);
	// Replace obsStateOfChildrenL by result times the lambda vector
	pL = resContraction.cwiseProduct(lambda);

	// Get scaling factors:
	double scalingProbL = x.getScalingFactorByEdgeId(iChildEdgeL);
	double scalingProbR = x.getScalingFactorByEdgeId(iChildEdgeR);
	double newScaling = scalingProbL+scalingProbR;
	newScaling += Likelihood::StateType::rescaleProbabilityVector(pL);
	// Set new scaling factor
	x.setScalingFactorByEdgeId(iChildEdgeL, newScaling);

	// Refresh the vec prob allocation: memory allocation of child L becomes memory allocation for parent edge
	x.updateVecProbForEdgeMapping(iChildEdgeL, iParentEdge);
	// Erase the unused vector: child R vector is now useless
	x.removeVecProbForEdge(iChildEdgeR);
}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::computeMassSpeciation(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {

	// Memorize who is a node event or who is a edge event and keep track of the involved vecProbs
	std::vector<int> vecProbToEdgeMapping(x.getVecProbToEdgeMapping());

	std::set<size_t> involvedVecProb;
	std::vector< std::pair<size_t, size_t> > involvedChildrenVecProbForNodes;
	for(size_t iN=0; iN<nodes.size(); iN++) { // For each node we find the two children edges
		assert(!nodes[iN]->isSampledAncestor() && !nodes[iN]->isExtinct() && !nodes[iN]->isExtant());

		std::vector<PS::Edge*> childrenEdges = nodes[iN]->getEdgesToChildren();
		std::vector<int>::iterator itFindL = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), childrenEdges[0]->getId());
		std::vector<int>::iterator itFindR = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), childrenEdges[1]->getId());
		assert(itFindL != vecProbToEdgeMapping.end() && itFindR != vecProbToEdgeMapping.end());

		// Iterator to indices
		size_t iVecProbL = std::distance(vecProbToEdgeMapping.begin(), itFindL);
		size_t iVecProbR = std::distance(vecProbToEdgeMapping.begin(), itFindR);

		// Keep track
		std::pair<size_t, size_t> edgesId = std::make_pair(iVecProbL, iVecProbR);
		involvedChildrenVecProbForNodes.push_back(edgesId);
		involvedVecProb.insert(iVecProbL);
		involvedVecProb.insert(iVecProbR);
	}

	// Get mass speciation event prob
	const Eigen::VectorXd& vecUpsilonOfT = ptrSynchEventContainer->getPtrMassSpeciation()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusUpsilonOfT = Eigen::VectorXd::Ones(vecUpsilonOfT.size()) - vecUpsilonOfT;

	// Get unobserved >REF< before update (carefull for openmp if all updated at once)
	Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();

	// First we compute for edges -  some work to do here for openmp implementation
	Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
	for(size_t iP=0; iP<x.size(); ++iP) {
		if(involvedVecProb.count(iP) > 0) continue; // We skip the vec prob involved in speciation events

		// Do the contractions
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iP);
		Eigen::VectorXd resContraction = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::MatrixXd >::ColXpr, Eigen::Ref< Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), p, u, t);

		// update p
		p = p.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContraction);

	}

	// Then we deal with the speciation events
	for(size_t iN=0; iN<nodes.size(); iN++) {
		std::pair<size_t, size_t> pair = involvedChildrenVecProbForNodes[iN];

		// Do the contraction
		Eigen::Ref< Eigen::MatrixXd >::ColXpr pL = obsStatesProb.col(pair.first);
		Eigen::Ref< Eigen::MatrixXd >::ColXpr pR = obsStatesProb.col(pair.second);

		// Do the contractions
		Eigen::VectorXd resContraction = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::MatrixXd >::ColXpr, Eigen::Ref< Eigen::MatrixXd >::ColXpr >(ptrTensorsContainer->getOmega(), pL, pR, t);


		// Storing in pLeft, pRight will have to be removed afterward
		pL = vecUpsilonOfT.cwiseProduct(resContraction);

	}

	// singleton:
	if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
		// Do the contractions
		Eigen::VectorXd resContractionS = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::VectorXd >, Eigen::Ref< Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), s, u, t);
		s = s.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContractionS);
	}

	// singleton no sampling:
	if(conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > sHat = x.getSingletonNoSamplingStateProb();
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();

		// Do the contractions
		Eigen::VectorXd resContractionSHat = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::VectorXd >, Eigen::Ref< Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), sHat, uHat, t);
		sHat = sHat.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContractionSHat);
	}

	// unobserved no sampling: must be updated last
	if(conditionalProbType == Conditions::ROOT_SURVIVAL || conditionalProbType == Conditions::ROOT_MRCA ||
	   conditionalProbType == Conditions::STEM_SURVIVAL || conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();

		// Do the contractions
		Eigen::VectorXd resContractionUHat = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::VectorXd >, Eigen::Ref< Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), uHat, uHat, t);
		uHat = vecOneMinusUpsilonOfT.cwiseProduct(uHat) + vecUpsilonOfT.cwiseProduct(resContractionUHat);
	}

	// Then we update the unobserved state probs: Must be done last
	// Do the contractions
	Eigen::VectorXd resContraction = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< Eigen::VectorXd >, Eigen::Ref< Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), u, u, t);
	u = vecOneMinusUpsilonOfT.cwiseProduct(u) + vecUpsilonOfT.cwiseProduct(resContraction);

	// Remove all the unused vec prob after speciation and update the edge / vec allocation
	for(size_t iN=0; iN<nodes.size(); iN++) {
		std::pair<size_t, size_t> pair = involvedChildrenVecProbForNodes[iN];
		x.removeVecProbForEdge(vecProbToEdgeMapping[pair.second]);
		x.updateVecProbForEdgeMapping(vecProbToEdgeMapping[pair.first], nodes[iN]->getEdgeToParent()->getId());
	}

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::computeMassExtinction(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {

	assert(nodes.empty());

	if(!ptrSynchEventContainer->getPtrMassExtinction()->areStateChangeInvolved()) { // No state change prob
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// First we compute for edges -  some work to do here for openmp implementation
		Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
		for(size_t iP=0; iP<x.size(); ++iP) {
			Eigen::Ref< Eigen::MatrixXd>::ColXpr p = obsStatesProb.col(iP) ;
			p = vecOneMinusGammaOfT.cwiseProduct(p);
		}

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
		u = vecGammaOfT +  vecOneMinusGammaOfT.cwiseProduct(u) ;

		// singleton:
		if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
			Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
			s = vecOneMinusGammaOfT.cwiseProduct(s);
		}

		// singleton no sampling:
		if(conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
			Eigen::Ref< Eigen::VectorXd > s_hat = x.getSingletonNoSamplingStateProb();
			s_hat = vecOneMinusGammaOfT.cwiseProduct(s_hat);
		}

		// unobserved no sampling:
		if(conditionalProbType == Conditions::ROOT_SURVIVAL || conditionalProbType == Conditions::ROOT_MRCA ||
		   conditionalProbType == Conditions::STEM_SURVIVAL || conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
			Eigen::Ref< Eigen::VectorXd > u_hat = x.getUnobservedNoSamplingStateProb();
			u_hat = vecGammaOfT + vecOneMinusGammaOfT.cwiseProduct(u_hat) ;
		}
	} else {
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		const Eigen::MatrixXd& matZetaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getStateChangeProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// First we compute for edges -  some work to do here for openmp implementation
		Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
		for(size_t iP=0; iP<x.size(); ++iP) {
			Eigen::Ref< Eigen::MatrixXd>::ColXpr p = obsStatesProb.col(iP) ;
			p = vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*p);
		}

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
		u = vecGammaOfT +  vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*u) ;

		// singleton:
		if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
			Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
			s = vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*s);
		}

		// singleton no sampling:
		if(conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
			Eigen::Ref< Eigen::VectorXd > s_hat = x.getSingletonNoSamplingStateProb();
			s_hat = vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*s_hat);
		}

		// unobserved no sampling:
		if(conditionalProbType == Conditions::ROOT_SURVIVAL || conditionalProbType == Conditions::ROOT_MRCA ||
		   conditionalProbType == Conditions::STEM_SURVIVAL || conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
			Eigen::Ref< Eigen::VectorXd > u_hat = x.getUnobservedNoSamplingStateProb();
			u_hat = vecGammaOfT + vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*u_hat) ;
		}
	}

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::computeMassSamplingEvent(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {

	// Define what to do on which node and with which vector
	std::vector<int> vecProbToEdgeMapping(x.getVecProbToEdgeMapping());

	// Keep track of what to do for each node/lineage
	std::vector<bool> isLineage(x.size(), true);
	std::vector< std::pair<size_t, size_t> > involvedExtinctionNodeAndVecProb;
	std::vector< std::pair<size_t, size_t> > involvedSampledAncNodeAndVecProb;
	for(size_t iN=0; iN<nodes.size(); iN++) { // For each node we defined what to do with it
		if(nodes[iN]->isExtinct()) {
			size_t iVecProb = x.allocateVecProbForEdge(nodes[iN]->getEdgeToParent()->getId());
			involvedExtinctionNodeAndVecProb.push_back(std::make_pair(iN, iVecProb));
			assert(iVecProb == isLineage.size());
			isLineage.push_back(false);
		} else if(nodes[iN]->isSampledAncestor()) {
			std::vector<int>::iterator itFind = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), nodes[iN]->getEdgesToChildren()[0]->getId());
			size_t iVecProb = std::distance(vecProbToEdgeMapping.begin(), itFind);
			involvedSampledAncNodeAndVecProb.push_back(std::make_pair(iN, iVecProb));
			isLineage[iVecProb] = false;
		} else {
			assert(false && "Node is neither extinct or sampled ancestor");
		}
	}

	// Unobserved probability
	Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();

	// Sampling probability
	const Eigen::VectorXd& vecRhoOfT = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusRhoOfT = Eigen::VectorXd::Ones(vecRhoOfT.size()) - vecRhoOfT;

	// Extinct nodes : init them
	Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
	for(size_t iN=0; iN<involvedExtinctionNodeAndVecProb.size(); ++iN) {
		size_t iNode = involvedExtinctionNodeAndVecProb[iN].first;
		size_t iVecProb = involvedExtinctionNodeAndVecProb[iN].second;
		PS::Node* extinctNode = nodes[iNode];
		assert(nodes[iNode]->isExtinct());

		// Get the observed data
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaId(extinctNode->getTaxaId());

		// Get the prob vector and initialize it
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iVecProb);
		p = probs.cwiseProduct(vecRhoOfT.cwiseProduct(u));
	}

	// Sampled ancestor node: update probability and vec prob allocation
	for(size_t iN=0; iN<involvedSampledAncNodeAndVecProb.size(); ++iN) {
		size_t iNode = involvedSampledAncNodeAndVecProb[iN].first;
		size_t iVecProb = involvedSampledAncNodeAndVecProb[iN].second;
		PS::Node* sampledAncestorNode = nodes[iNode];
		assert(nodes[iNode]->isSampledAncestor());


		// keep track of edges id
		size_t iBeforeEdge = sampledAncestorNode->getEdgesToChildren().front()->getId();
		size_t iAfterEdge = sampledAncestorNode->getEdgeToParent()->getId();

		// Get the observed data
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaId(sampledAncestorNode->getTaxaId());

		// Get the vector as a function of the before edge ID
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iVecProb);
		Eigen::VectorXd tmpP = p; // Copy
		p = probs.cwiseProduct(vecRhoOfT.cwiseProduct(tmpP));

		// Update the Edge to prob vector mapping
		x.updateVecProbForEdgeMapping(iBeforeEdge, iAfterEdge);
	}

	// Lineage case:
	assert(isLineage.size() == x.size());
	for(size_t iP=0; iP<x.size(); ++iP) {
		if(!isLineage[iP]) continue; // Skip already treated nodes
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iP);
		p = vecOneMinusRhoOfT.cwiseProduct(p);
	}

	// Singleton case:
	if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
		s = vecOneMinusRhoOfT.cwiseProduct(s) + vecRhoOfT.cwiseProduct(u);
	}

	// Unobserved
	u = vecOneMinusRhoOfT.cwiseProduct(u);

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::computeMassDestrSamplingEvent(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {

	// Define what to do on which node and with which vector
	std::vector<int> vecProbToEdgeMapping(x.getVecProbToEdgeMapping());

	// Keep track of what to do for each node/lineage
	std::vector<bool> isLineage(x.size(), true);
	std::vector< std::pair<size_t, size_t> > involvedExtinctionNodeAndVecProb;
	for(size_t iN=0; iN<nodes.size(); iN++) { // For each node we defined what to do with it
		if(nodes[iN]->isExtinct()) {
			size_t iVecProb = x.allocateVecProbForEdge(nodes[iN]->getEdgeToParent()->getId());
			involvedExtinctionNodeAndVecProb.push_back(std::make_pair(iN, iVecProb));
			assert(iVecProb == isLineage.size());
			isLineage.push_back(false);
		}  else {
			assert(false && "Node is not extinct");
		}
	}

	// Sampling probability
	const Eigen::VectorXd& vecXiOfT = ptrSynchEventContainer->getPtrMassDestrSampling()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusXiOfT = Eigen::VectorXd::Ones(vecXiOfT.size()) - vecXiOfT;

	// Extinct nodes : init them
	Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
	for(size_t iN=0; iN<involvedExtinctionNodeAndVecProb.size(); ++iN) {
		size_t iNode = involvedExtinctionNodeAndVecProb[iN].first;
		size_t iVecProb = involvedExtinctionNodeAndVecProb[iN].second;
		PS::Node* extinctNode = nodes[iNode];
		assert(nodes[iNode]->isExtinct());

		// Get the observed data
		const Eigen::VectorXd &probs = ptrData->getProbForTaxaId(extinctNode->getTaxaId());

		// Get the prob vector and initialize it
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iVecProb);
		p = vecXiOfT.cwiseProduct(probs);
	}

	// Lineage case:
	assert(isLineage.size() == x.size());
	for(size_t iP=0; iP<x.size(); ++iP) {
		if(!isLineage[iP]) continue; // Skip already treated nodes
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iP);
		p = vecOneMinusXiOfT.cwiseProduct(p);
	}

	// Singleton case:
	if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
		s = vecOneMinusXiOfT.cwiseProduct(s) + vecXiOfT;;
	}


	// Unobserved probability
	Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
	u = vecOneMinusXiOfT.cwiseProduct(u);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
void EigenKernels<conditionalProbType, withCladoEvents>::rescaleRequestedBranchesProbabilities(const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {

	for(size_t iN=0; iN<nodes.size(); ++iN) {

		// Get children vec prob
		double iEdge = nodes[iN]->getEdgeToParent()->getId();
		Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
		size_t iP = x.findVecProbIdForEdgeId(iEdge);
		Eigen::Ref< Eigen::MatrixXd >::ColXpr p = obsStatesProb.col(iP);
		double scalingFactor = x.getScalingFactorByEdgeId(iEdge);
		scalingFactor += Likelihood::StateType::rescaleProbabilityVector(p);

		// Set new scaling factor
		x.setScalingFactorByEdgeId(iEdge, scalingFactor);
	}

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
double EigenKernels<conditionalProbType, withCladoEvents>::computeLogLikelihood(double t, PS::Node* originNode, Eigen::VectorXd &prior, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x) {
	assert(originNode->isOriginNode());
	assert(x.size() == 1);

	Eigen::Ref< Eigen::MatrixXd > obsStatesProb = x.getObservedStateProb();
	Eigen::Ref< Eigen::MatrixXd >::ColXpr p= obsStatesProb.col(0);
	Eigen::Ref< Eigen::VectorXd >u = x.getUnobservedStateProb();
	assert(p.size() == prior.size());

	Eigen::VectorXd s = Eigen::VectorXd::Ones(prior.size()) - u;
	double unscaledPData = prior.dot(p);
	double logPData = log(unscaledPData) + x.getScalingFactorByVecPos(0);
	double pSurvival = 1.0;

	if(conditionalProbType == Likelihood::Conditions::TIME ) {
		// Do nothing
	} else if (conditionalProbType == Likelihood::Conditions::ROOT_SURVIVAL ) {
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();
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
	} else if (conditionalProbType == Likelihood::Conditions::ROOT_MRCA ) {
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();
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
	} else if (conditionalProbType == Likelihood::Conditions::ROOT_SAMPLING ) {
		Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
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
	} else if (conditionalProbType == Likelihood::Conditions::STEM_SURVIVAL ) {
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();
		Eigen::VectorXd vecOneMinusUHat = Eigen::VectorXd::Ones(uHat.size()) - uHat;
		pSurvival = prior.dot(vecOneMinusUHat);
	} else if (conditionalProbType == Likelihood::Conditions::STEM_ONE_SAMPLE ) {
		Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
		Eigen::VectorXd vecOneMinusU = Eigen::VectorXd::Ones(u.size()) - u;
		pSurvival = prior.dot(vecOneMinusU);
	} else if (conditionalProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES ) {
		Eigen::Ref< Eigen::VectorXd > uHat = x.getUnobservedNoSamplingStateProb();
		Eigen::Ref< Eigen::VectorXd > sHat = x.getSingletonNoSamplingStateProb();
		Eigen::VectorXd vecOneMinusUHatSHat = Eigen::VectorXd::Ones(uHat.size()) - uHat ;
		vecOneMinusUHatSHat -= sHat;
		pSurvival = prior.dot(vecOneMinusUHatSHat);
	} else if (conditionalProbType == Likelihood::Conditions::STEM_TWO_SAMPLES ) {
		Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
		Eigen::Ref< Eigen::VectorXd > s = x.getSingletonStateProb();
		Eigen::VectorXd vecOneMinusUS = Eigen::VectorXd::Ones(u.size()) - u;
		vecOneMinusUS -= s;
		pSurvival = prior.dot(vecOneMinusUS);
	} else if (conditionalProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
		Eigen::Ref< Eigen::VectorXd > u = x.getUnobservedStateProb();
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

	return logPData - log(pSurvival);
}


} /* namespace OpenMP */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif //defined(_OPENMP)
#endif // LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENKERNELS_DEF_H_
