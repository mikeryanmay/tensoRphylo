/*
 * AsynchronousApproximator.cpp
 *
 *  Created on: Apr 17, 2020
 *      Author: meyerx
 */

#include "AsynchronousApproximator.h"

#include "Data/Structure/Tree.h"
#include "IncLikelihoodApproximator.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Utils/Output/OutputManager.h"
#include "Likelihood/Scheduler/DAG/IncDAG.h"

namespace Likelihood {
namespace Approximator {

AsynchronousApproximator::AsynchronousApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
		 Conditions::conditionalProbability_t aConditionType,
		 Phylogeny::Data::ContainerSharedPtr aPtrData,
		 Scheduler::SchedulerSharedPtr aPtrScheduler,
		 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		 Tensor::ContainerSharedPtr aPtrTensorCont) :
						BaseApproximator(aIntScheme, aConditionType, aPtrData,
								aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
						ptrSApproxManager(new Specialized::SpecializedApproxManager(aConditionType, ptrScheduler, ptrSyncEventsCont, ptrTensorCont)),
						ptrDAG(new Likelihood::Scheduler::DAG::DAG(ptrScheduler)),
						ptrSchedDAG(new Likelihood::Scheduler::DAG::SchedulerDAG(ptrDAG)) {

}

AsynchronousApproximator::~AsynchronousApproximator() {

}

double AsynchronousApproximator::approximateLogLikelihood() {

	logLikelihood = -std::numeric_limits<double>::max();

	// Enable specialized implementation operations: e.g., init state vectors
	doPreProcessingSteps();

	// Check the event validity and return 0 if an event is impossible
	if(!areEventsPossible()) {
		std::cerr << "TensorPhylo detected invalid combination of parameters. Hint: simultaneous events that are incompatible are defined (logLik=-inf)." << std::endl;
		return logLikelihood;
	}

	// First init specialized
	ptrSApproxManager->compute();

	{ // Then use DAG

		using namespace Likelihood::Scheduler::DAG;

		while(!ptrSchedDAG->isDone()) {
			NodeDAG* task = ptrSchedDAG->getNextTask();
			assert(task != NULL);
			assert(task->isReady());
			//std::cout << task->toString() << std::endl;

			if(task->isEventTask()) {
				doEventStep(task);
			} else if(task->isIntegrationTask()) {
				doIntegrationStep(task);
			} else {
				assert(false && "Unknown type of task.");
			}
			ptrSchedDAG->signalDone(task);
		}

	}

	// Enable specialized implementation operations: e.g., cleanup
	doPostProcessingSteps();

	// We should be good
	if(applyTreeCorrection) {
		double logLikelihoodCorrection = ptrScheduler->getPtrTree()->getLogLikelihoodCorrection();
		logLikelihood += logLikelihoodCorrection;
	}
	return logLikelihood;

}

} /* namespace Approximator */
} /* namespace Likelihood */
