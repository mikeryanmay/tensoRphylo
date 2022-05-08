/*
 * EigenKernelsPairedUS.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "EigenKernelsPairedUS.h"

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

EigenKernelsPairedUS::EigenKernelsPairedUS(Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
							SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
							Tensor::ContainerSharedPtr aPtrTensorsContainer) :
		PROCESS_TYPE(aExtantProcessType),
		ptrSynchEventContainer(aPtrSynchEventContainer),
		ptrTensorsContainer(aPtrTensorsContainer) {
}

EigenKernelsPairedUS::~EigenKernelsPairedUS() {
}


void EigenKernelsPairedUS::setInitialCondition(Likelihood::StateType::Matrix::EigenState &x) {

	x.resize(2);

	// Get initial mass sampling event
	const Eigen::VectorXd& vecRhoOf0 = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(0.);

	// Init the unobserved probs
	Eigen::VectorXd& unobservedStatesProb = x.getStateProb(0);
	unobservedStatesProb.setOnes();
	unobservedStatesProb -= vecRhoOf0;

	// Init the unobserved probs
	Eigen::VectorXd& s = x.getStateProb(1);
	s = vecRhoOf0;

}

void EigenKernelsPairedUS::computeMassSpeciation(double t, Likelihood::StateType::Matrix::EigenState &x) {

	const Eigen::VectorXd& vecUpsilonOfT = ptrSynchEventContainer->getPtrMassSpeciation()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusUpsilonOfT = Eigen::VectorXd::Ones(vecUpsilonOfT.size()) - vecUpsilonOfT;

	// Get U or u hat
	Eigen::VectorXd& u = x.getStateProb(0);

	Eigen::VectorXd& s = x.getStateProb(1);
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

	// Then we update the unobserved state probs: Must be done last
	if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
		Eigen::VectorXd resContraction = u.cwiseProduct(u);
		u = vecOneMinusUpsilonOfT.cwiseProduct(u) + vecUpsilonOfT.cwiseProduct(resContraction);
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		// Do the contraction
		Eigen::VectorXd resContraction = computeTensorContractionVector(omega, u, u);
		u = vecOneMinusUpsilonOfT.cwiseProduct(u) + vecUpsilonOfT.cwiseProduct(resContraction);
	} else {
		Tensor::tensor_t &omega = ptrTensorsContainer->getOmega()->getTensor(t);
		// Do the contraction
		Eigen::VectorXd resContraction = computeTensorContractionVector(omega, u, u);
		u = vecOneMinusUpsilonOfT.cwiseProduct(u) + vecUpsilonOfT.cwiseProduct(resContraction);
	}


}

void EigenKernelsPairedUS::computeMassExtinction(double t,  Likelihood::StateType::Matrix::EigenState &x) {

	if(!ptrSynchEventContainer->getPtrMassExtinction()->areStateChangeInvolved()) { // No state change prob
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::VectorXd& u = x.getStateProb(0);
		u = vecGammaOfT +  vecOneMinusGammaOfT.cwiseProduct(u) ;
		// Singleton
		Eigen::VectorXd& s = x.getStateProb(1);
		s = vecOneMinusGammaOfT.cwiseProduct(s);

	} else {
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		const Eigen::MatrixXd& matZetaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getStateChangeProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::VectorXd& u = x.getStateProb(0);
		u = vecGammaOfT +  vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*u) ;
		// Singleton
		Eigen::VectorXd& s = x.getStateProb(1);
		s = vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*s);

	}

}

void EigenKernelsPairedUS::computeMassSamplingEvent(double t, Likelihood::StateType::Matrix::EigenState &x) {

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {

		// Unobserved probability
		Eigen::VectorXd& u = x.getStateProb(0);

		// Sampling probability
		const Eigen::VectorXd& vecRhoOfT = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusRhoOfT = Eigen::VectorXd::Ones(vecRhoOfT.size()) - vecRhoOfT;

		// singleton probability
		Eigen::VectorXd& s = x.getStateProb(1);

		// singleton
		s = vecOneMinusRhoOfT.cwiseProduct(s) + vecRhoOfT.cwiseProduct(u);

		// Unobserved
		u = vecOneMinusRhoOfT.cwiseProduct(u);
	}

}

void EigenKernelsPairedUS::computeMassDestrSamplingEvent(double t, Likelihood::StateType::Matrix::EigenState &x) {

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		// Sampling probability
		const Eigen::VectorXd& vecXiOfT = ptrSynchEventContainer->getPtrMassDestrSampling()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusXiOfT = Eigen::VectorXd::Ones(vecXiOfT.size()) - vecXiOfT;

		// Unobserved probability
		Eigen::VectorXd& u = x.getStateProb(0);
		u = vecOneMinusXiOfT.cwiseProduct(u);

		// singleton probability
		Eigen::VectorXd& s = x.getStateProb(1);
		s = vecOneMinusXiOfT.cwiseProduct(s) + vecXiOfT;

	}
}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
