/*
 * EigenKernelsUnobserved.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "EigenKernelsUnobserved.h"

#include "../EigenUtils.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/StateTypes/Utils.h"
#include "Likelihood/Scheduler/IncScheduler.h"

#include <Eigen/Core>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

EigenKernelsUnobserved::EigenKernelsUnobserved(
							Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
							SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
							Tensor::ContainerSharedPtr aPtrTensorsContainer) :
									PROCESS_TYPE(aExtantProcessType),
									ptrSynchEventContainer(aPtrSynchEventContainer),
									ptrTensorsContainer(aPtrTensorsContainer) {
}

EigenKernelsUnobserved::~EigenKernelsUnobserved() {
}


void EigenKernelsUnobserved::setInitialCondition(Likelihood::StateType::Vector::EigenState &x) {

	// Get initial mass sampling event
	const Eigen::VectorXd& vecRhoOf0 = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(0.);

	// Init the unobserved probs
	Eigen::VectorXd& unobservedStatesProb = x.getStateProb();
	unobservedStatesProb.setOnes();
	unobservedStatesProb -= vecRhoOf0;

}

void EigenKernelsUnobserved::computeMassSpeciation(double t, Likelihood::StateType::Vector::EigenState &x) {

	Eigen::VectorXd& u = x.getStateProb();

	const Eigen::VectorXd& vecUpsilonOfT = ptrSynchEventContainer->getPtrMassSpeciation()->getEventProbability(t);
	Eigen::VectorXd vecOneMinusUpsilonOfT = Eigen::VectorXd::Ones(vecUpsilonOfT.size()) - vecUpsilonOfT;

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

void EigenKernelsUnobserved::computeMassExtinction(double t,  Likelihood::StateType::Vector::EigenState &x) {

	if(!ptrSynchEventContainer->getPtrMassExtinction()->areStateChangeInvolved()) { // No state change prob
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::VectorXd& u = x.getStateProb();
		u = vecGammaOfT +  vecOneMinusGammaOfT.cwiseProduct(u) ;

	} else {
		// Get mass extinction event probs
		const Eigen::VectorXd& vecGammaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getEventProbability(t);
		const Eigen::MatrixXd& matZetaOfT = ptrSynchEventContainer->getPtrMassExtinction()->getStateChangeProbability(t);
		Eigen::VectorXd vecOneMinusGammaOfT = Eigen::VectorXd::Ones(vecGammaOfT.size()) - vecGammaOfT;

		// Get unobserved >REF< before update (carefull for openmp if all updated at once)
		Eigen::VectorXd& u = x.getStateProb();
		u = vecGammaOfT +  vecOneMinusGammaOfT.cwiseProduct(matZetaOfT*u) ;

	}

}

void EigenKernelsUnobserved::computeMassSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x) {

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		// Unobserved probability
		Eigen::VectorXd& u = x.getStateProb();

		// Sampling probability
		const Eigen::VectorXd& vecRhoOfT = ptrSynchEventContainer->getPtrMassSampling()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusRhoOfT = Eigen::VectorXd::Ones(vecRhoOfT.size()) - vecRhoOfT;

		// Unobserved
		u = vecOneMinusRhoOfT.cwiseProduct(u);
	}

}

void EigenKernelsUnobserved::computeMassDestrSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x) {

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		// Sampling probability
		const Eigen::VectorXd& vecXiOfT = ptrSynchEventContainer->getPtrMassDestrSampling()->getEventProbability(t);
		Eigen::VectorXd vecOneMinusXiOfT = Eigen::VectorXd::Ones(vecXiOfT.size()) - vecXiOfT;

		// Unobserved probability
		Eigen::VectorXd& u = x.getStateProb();
		u = vecOneMinusXiOfT.cwiseProduct(u);

	}
}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
