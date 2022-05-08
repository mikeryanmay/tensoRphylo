/*
 * SynchronousApproximator.cpp
 *
 *  Created on: Apr 17, 2020
 *      Author: meyerx
 */

#include "SynchronousApproximator.h"

#include "Data/Structure/Tree.h"
#include "Likelihood/Scheduler/IncScheduler.h"

namespace Likelihood {
namespace Approximator {

SynchronousApproximator::SynchronousApproximator(
		Likelihood::Integrator::integrationScheme_t aIntScheme,
		Conditions::conditionalProbability_t aConditionType,
		Phylogeny::Data::ContainerSharedPtr aPtrData,
		Scheduler::SchedulerSharedPtr aPtrScheduler,
		SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		Tensor::ContainerSharedPtr aPtrTensorCont) :
				BaseApproximator(aIntScheme, aConditionType, aPtrData,
						aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont) {
}

SynchronousApproximator::~SynchronousApproximator() {
}

double SynchronousApproximator::approximateLogLikelihood() {

	// By now we don't care if it has been updated
	ptrScheduler->clearHasBeenUpdatedFlag();

	logLikelihood = -std::numeric_limits<double>::max();

	// Enable specialized implementation operations: e.g., init state vectors
	doPreProcessingSteps();

	// Check the event validity and return 0 if an event is impossible
	if(!areEventsPossible()) {
		std::cerr << "TensorPhylo detected invalid combination of parameters. Hint: simultaneous events that are incompatible are defined (logLik=-inf)." << std::endl;
		return logLikelihood;
	}

	size_t iEvent = 0;
	// Loop over integration/event steps until we reach the final event
	while(!ptrScheduler->getEvents()[iEvent]->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT)) {
		// Do event related computations
		//std::cout << ptrScheduler->getEvents()[iEvent]->toString() << std::endl;
		doEventStep(iEvent);
		// Integrate over the edges
		doIntegrationStep(iEvent);
		// Next event
		iEvent ++;
	}

	// Final event compute the likelihood
	doEventStep(iEvent);

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
