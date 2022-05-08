/*
 * AdaptiveUnobserved.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "AdaptiveUnobserved.h"

#include <cmath>

#include "Tensor/Container.h"
#include "Data/Structure/Node.h"
#include "../BaseApproximator.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.h"


namespace Likelihood {
namespace Approximator {
namespace Specialized {

AdaptiveUnobserved::AdaptiveUnobserved(
		processType_t aExtantProcessType,
		Scheduler::SchedulerSharedPtr aPtrScheduler,
		SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		Tensor::ContainerSharedPtr aPtrTensorCont) :
	BaseSpecialized(aExtantProcessType, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
	kernels(BaseSpecialized::PROCESS_TYPE, aPtrSyncEventsCont, aPtrTensorCont),
	intKernel(BaseSpecialized::PROCESS_TYPE, aPtrTensorCont) {
}

AdaptiveUnobserved::~AdaptiveUnobserved() {
}

void AdaptiveUnobserved::doIntegrationStep(double startTime, double endTime) {
	// We can get the layer of edges using the scheduler
	// However, it's not even needed as the kernels already deal with that
	// After a doEventStep, the state has already been resized to fit the next set of edges

	if(endTime == startTime) return;

	if(endTime != ptrScheduler->getEvents().back()->getTime()) {
		endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	}

	Likelihood::Integrator::Base<StateType, IntegrationKernelsType, OperationType>* ptrIntegrator;
	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<StateType, IntegrationKernelsType, OperationType >(BaseApproximator::DEFAULT_ABS_TOLERANCE,
			BaseApproximator::DEFAULT_REL_TOLERANCE,
			BaseApproximator::DEFAULT_DELTA_T,
			Integrator::RUNGE_KUTTA_DOPRI5);

	ptrIntegrator->integrate(startTime, endTime, probState, intKernel);

	delete ptrIntegrator;

}

void AdaptiveUnobserved::doEventStep(size_t iEvent) {

	Likelihood::Scheduler::Event* event = ptrScheduler->getEvents()[iEvent];

	if(event->checkEvent(Likelihood::Scheduler::PRESENT_TIME_EVENT)) {
		kernels.setInitialCondition(probState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT)) {
		kernels.computeMassSpeciation(event->getTime(), probState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT)) {
		kernels.computeMassExtinction(event->getTime(), probState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT)) {
		kernels.computeMassSamplingEvent(event->getTime(), probState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		kernels.computeMassDestrSamplingEvent(event->getTime(), probState);
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
