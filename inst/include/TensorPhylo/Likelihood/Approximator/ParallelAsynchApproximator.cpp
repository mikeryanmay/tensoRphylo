/*
 * ParallelAsynchApproximator.cpp
 *
 *  Created on: Apr 17, 2020
 *      Author: meyerx
 */

#include "ParallelAsynchApproximator.h"

#include "Utils/Parallel/Manager.h"
#if defined(_OPENMP)

#include "Data/Structure/Tree.h"
#include "IncLikelihoodApproximator.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Utils/Output/OutputManager.h"
#include "Likelihood/Scheduler/DAG/IncDAG.h"

namespace Likelihood {
namespace Approximator {

ParallelAsynchApproximator::ParallelAsynchApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
		 Conditions::conditionalProbability_t aConditionType,
		 Phylogeny::Data::ContainerSharedPtr aPtrData,
		 Scheduler::SchedulerSharedPtr aPtrScheduler,
		 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		 Tensor::ContainerSharedPtr aPtrTensorCont) :
						BaseApproximator(aIntScheme, aConditionType, aPtrData,
								aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
						ptrSApproxManager(new Specialized::SpecializedApproxManager(aConditionType, ptrScheduler, ptrSyncEventsCont, ptrTensorCont)),
						ptrDAG(new Likelihood::Scheduler::DAG::DAG(ptrScheduler)),
						ptrSchedDAG(new Likelihood::Scheduler::DAG::ParallelSchedulerDAG(Utils::Parallel::Manager::getInstance()->getNThread(), ptrDAG)) {

}

ParallelAsynchApproximator::~ParallelAsynchApproximator() {

}

double ParallelAsynchApproximator::approximateLogLikelihood() {

	logLikelihood = -std::numeric_limits<double>::max();
	if(ptrSchedDAG->getNThread() != Utils::Parallel::Manager::getInstance()->getNThread()) {
		ptrSchedDAG.reset(new Likelihood::Scheduler::DAG::ParallelSchedulerDAG(Utils::Parallel::Manager::getInstance()->getNThread(), ptrDAG));
	}

	// Enable specialized implementation operations: e.g., init state vectors
	doPreProcessingSteps();

	// Check the event validity and return 0 if an event is impossible
	if(!areEventsPossible()) {
		std::cerr << "TensorPhylo detected invalid combination of parameters. Hint: simultaneous events that are incompatible are defined (logLik=-inf)." << std::endl;
		return logLikelihood;
	}


	// First init specialized
	ptrSApproxManager->compute();


	#pragma omp parallel
	{ // Then use DAG

		#pragma omp task
		{
			using namespace Likelihood::Scheduler::DAG;

			size_t iThread = omp_get_thread_num();

			NodeDAG* task = ptrSchedDAG->getNextTask(iThread);
			while(!ptrSchedDAG->isDone() && task != NULL) {
				assert(task != NULL);
				assert(task->isReady());

				if(task->isEventTask()) {
					doEventStep(task);
				} else if(task->isIntegrationTask()) {
					doIntegrationStep(task);
				} else {
					assert(false && "Unknown type of task.");
				}

				// Trying to get a task without having to go through a critical section
				bool success = ptrSchedDAG->tryNonBlockingSignalAndPop(iThread, task);
				if(!success) { // failure, trying through a critical section
					task = ptrSchedDAG->getNextTask(iThread);
				}
			}
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

#endif
