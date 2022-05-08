/*
 * BaseSpecialized.cpp
 *
 *  Created on: Apr 16, 2020
 *      Author: meyerx
 */

#include "BaseSpecialized.h"

#include <Eigen/Core>

#include "Tensor/Container.h"
#include "Utils/MemoryPool/EigenCPU.h"
#include "Likelihood/Scheduler/Event.h"
#include "Likelihood/Scheduler/BaseScheduler.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

BaseSpecialized::BaseSpecialized(
		processType_t aExtantProcessType,
		Scheduler::SchedulerSharedPtr aPtrScheduler,
		SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		Tensor::ContainerSharedPtr aPtrTensorCont) :
		PROCESS_TYPE(aExtantProcessType),
	ptrScheduler(aPtrScheduler) {

	ready = false;
	Utils::MemoryPool::eigenCPU().setNCategories(aPtrTensorCont->getNumberOfState());
}

BaseSpecialized::~BaseSpecialized() {
}

const Eigen::VectorXd& BaseSpecialized::getFinalProbability() const {
	assert(ready && "compute() should be called first.");
	return probState.getStateProb();
}

bool BaseSpecialized::isExtantProcess() const {
	return PROCESS_TYPE == EXTANT_PROCESS;
}

bool BaseSpecialized::isFullProcess() const {
	return PROCESS_TYPE == FULL_PROCESS;
}

processType_t BaseSpecialized::getProcessType() const {
	return PROCESS_TYPE;
}

void BaseSpecialized::compute() {

	// Enable specialized implementation operations: e.g., init state vectors
	doPreprocessingStep();

	size_t iEvent = 0;
	bool done = false;
	while(!done) {
		// Do event related computations
		//std::cout << "Current event : " << ptrScheduler->getEvents()[iEvent]->toString() << std::endl;
		doEventStep(iEvent);
		double curEventTime = ptrScheduler->getEvents()[iEvent]->getTime();

		// Advance until next valid event (for unobserved)
		bool foundNext = false;
		while(!foundNext) {
			iEvent++;
			//std::cout << "Trying next event : " << ptrScheduler->getEvents()[iEvent]->toString() << std::endl;
			foundNext = ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT) ||
					 ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT) ||
				     ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT) ||
					 ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT) ||
					 ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT) ||
					 ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RATE_SHIFT);
		}
		//std::cout << "Next event : " << ptrScheduler->getEvents()[iEvent]->toString() << std::endl;
		//std::cout << "--------------------------------------------------------" << std::endl;
		double nextEventTime = ptrScheduler->getEvents()[iEvent]->getTime();

		// Integrate until next event
		doIntegrationStep(curEventTime, nextEventTime);

		// Stopping condition
		done = ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT);
	}

	ready = true;
}

void BaseSpecialized::doPreprocessingStep() {
	probState = StateType(); // reset state
}

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */
