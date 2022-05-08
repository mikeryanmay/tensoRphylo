/*
 * AdaptivePairedUS.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "AdaptivePairedUS.h"

#include <cmath>

#include "Data/Structure/Node.h"
#include "Tensor/Container.h"
#include "../BaseApproximator.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.h"


namespace Likelihood {
namespace Approximator {
namespace Specialized {

AdaptivePairedUS::AdaptivePairedUS(processType_t aExtantProcessType,
		Scheduler::SchedulerSharedPtr aPtrScheduler,
		SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		Tensor::ContainerSharedPtr aPtrTensorCont) :
			BaseSpecialized(aExtantProcessType, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
			kernels(BaseSpecialized::PROCESS_TYPE, aPtrSyncEventsCont, aPtrTensorCont),
			intKernel(BaseSpecialized::PROCESS_TYPE, aPtrTensorCont) {

}

AdaptivePairedUS::~AdaptivePairedUS() {
}

const Eigen::VectorXd& AdaptivePairedUS::getU() const {
	assert(matrixProbState.size() == 2);
	return matrixProbState.getStateProb(0);
}

const Eigen::VectorXd& AdaptivePairedUS::getS() const {
	assert(matrixProbState.size() == 2);
	return matrixProbState.getStateProb(1);
}

void AdaptivePairedUS::doPreprocessingStep() {
	matrixProbState = StateType();
}

void AdaptivePairedUS::doIntegrationStep(double startTime, double endTime) {
	// We can get the layer of edges using the scheduler
	// However, it's not even needed as the kernels already deal with that
	// After a doEventStep, the state has already been resized to fit the next set of edges

	if(endTime == startTime) return;

	if(endTime != ptrScheduler->getEvents().back()->getTime()) {
		endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	}

	Likelihood::Integrator::Base<StateType, IntegrationKernelType, OperationType>* ptrIntegrator;
	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<StateType, IntegrationKernelType, OperationType >(BaseApproximator::DEFAULT_ABS_TOLERANCE,
			BaseApproximator::DEFAULT_REL_TOLERANCE,
			BaseApproximator::DEFAULT_DELTA_T,
			Integrator::RUNGE_KUTTA_DOPRI5);

	ptrIntegrator->integrate(startTime, endTime, matrixProbState, intKernel);

	delete ptrIntegrator;

}

void AdaptivePairedUS::doEventStep(size_t iEvent) {

	Likelihood::Scheduler::Event* event = ptrScheduler->getEvents()[iEvent];

	if(event->checkEvent(Likelihood::Scheduler::PRESENT_TIME_EVENT)) {
		kernels.setInitialCondition(matrixProbState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT)) {
		kernels.computeMassSpeciation(event->getTime(), matrixProbState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT)) {
		kernels.computeMassExtinction(event->getTime(), matrixProbState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT)) {
		kernels.computeMassSamplingEvent(event->getTime(), matrixProbState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		kernels.computeMassDestrSamplingEvent(event->getTime(), matrixProbState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RATE_SHIFT)) {
		// DO NOTHING BUT ENSURE THAT THE NUMBERICAL INTEGRATOR DEAL ONLY WITH CONTIUNOUS FUNCTIONS
	} else {
		std::cout << event->toString() << std::endl;
		std::cout << event->getNodes().front()->toString() << std::endl;
		assert(false && "Event is not implemented.");
	}

}

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */
