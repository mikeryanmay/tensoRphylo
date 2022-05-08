/*
 * EigenKernelsSingleton.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "EigenKernelsSingleton.h"

#include "../EigenUtils.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/StateTypes/Utils.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/Approximator/Specialized/DenseUnobserved.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"

#include <Eigen/Core>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

EigenKernelsSingleton::EigenKernelsSingleton(Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
							Likelihood::Approximator::Specialized::DenseUnobservedSharedPtr aPtrUnobsProb,
							SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
							Tensor::ContainerSharedPtr aPtrTensorsContainer) :
		PROCESS_TYPE(aExtantProcessType),
		ptrDenseUnobs(aPtrUnobsProb),
		ptrSynchEventContainer(aPtrSynchEventContainer),
		ptrTensorsContainer(aPtrTensorsContainer) {
	assert(PROCESS_TYPE == ptrDenseUnobs->getProcessType());
}

EigenKernelsSingleton::~EigenKernelsSingleton() {
}


void EigenKernelsSingleton::setInitialCondition(Likelihood::StateType::Vector::EigenState &x) {

	// Get initial mass sampling event
	const Eigen::VectorXd& vecRhoOf0 = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(0.);

	// Init the unobserved probs
	Eigen::VectorXd& s = x.getStateProb();
	s = vecRhoOf0;

}

void EigenKernelsSingleton::computeMassSpeciation(double t, Likelihood::StateType::Vector::EigenState &x) {

	const Eigen::VectorXd& vecUpsilonOfT = ptrSynchEventContainer->getPtrMassSpeciation()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusUpsilonOfT = Eigen::VectorXd::Ones(vecUpsilonOfT.size()) - vecUpsilonOfT;

	// Get U or u hat
	Eigen::VectorXd u;
	ptrDenseUnobs->defineProbabilities(t, u);

	Eigen::VectorXd& s = x.getStateProb();
	if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
		Eigen::VectorXd resContraction = s.cwiseProduct(u);
		s = s.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContraction);
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		// Do the contractions
		Eigen::VectorXd resContractionS = computeTensorContractionVector(omega, s, u);
		s = s.cwiseProduct(vecOneMinusUpsilonOfT) + 2.*vecUpsilonOfT.cwiseProduct(resContractionS);
	} else {
		Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
		// Do the contractions
		Eigen::VectorXd resContractionS = computeTensorContractionVector(omega, s, u);
		s = s.cwiseProduct(vecOneMinusUpsilonOfT) + 2.* vecUpsilonOfT.cwiseProduct(resContractionS);
	}

}

void EigenKernelsSingleton::computeMassExtinction(double t,  Likelihood::StateType::Vector::EigenState &x) {

	if(!ptrSynchEventContainer->getPtrMassExtinction()->areStateChangeInvolved()) { // No state change prob
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::VectorXd& s = x.getStateProb();
		s = vecOneMinusGammaOfT.cwiseProduct(s);
	} else {
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		const Eigen::MatrixXd& matZetaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getStateChangeProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::VectorXd& s = x.getStateProb();
		s = vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*s);
	}

}

void EigenKernelsSingleton::computeMassSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x) {

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		// Unobserved probability
		Eigen::VectorXd& s = x.getStateProb();

		// Sampling probability
		const Eigen::VectorXd& vecRhoOfT = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusRhoOfT = Eigen::VectorXd::Ones(vecRhoOfT.size()) - vecRhoOfT;

		// Compute u
		Eigen::VectorXd u;
		ptrDenseUnobs->defineProbabilities(t, u);

		// Unobserved
		s = vecOneMinusRhoOfT.cwiseProduct(s) + vecRhoOfT.cwiseProduct(u);
	}

}

void EigenKernelsSingleton::computeMassDestrSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x) {

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		// Sampling probability
		const Eigen::VectorXd& vecXiOfT = ptrSynchEventContainer->getPtrMassDestrSampling()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusXiOfT = Eigen::VectorXd::Ones(vecXiOfT.size()) - vecXiOfT;

		// Unobserved probability
		Eigen::VectorXd& s = x.getStateProb();
		s = vecOneMinusXiOfT.cwiseProduct(s) + vecXiOfT;

	}
}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
