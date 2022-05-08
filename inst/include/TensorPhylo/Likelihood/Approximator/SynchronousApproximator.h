/*
 * SynchronousApproximator.h
 *
 *  Created on: Apr 17, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_SYNCHRONOUSAPPROXIMATOR_H_
#define LIKELIHOOD_SYNCHRONOUSAPPROXIMATOR_H_

#include "BaseApproximator.h"

namespace Likelihood {
namespace Approximator {

class SynchronousApproximator: public BaseApproximator {
public:
	SynchronousApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
			 Conditions::conditionalProbability_t aConditionType,
			 Phylogeny::Data::ContainerSharedPtr aPtrData,
			 Scheduler::SchedulerSharedPtr aPtrScheduler,
			 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			 Tensor::ContainerSharedPtr aPtrTensorCont);
	virtual ~SynchronousApproximator();


	double approximateLogLikelihood();

protected:

	virtual void doPreProcessingSteps() = 0;
	virtual void doIntegrationStep(size_t iEdgesLayer) = 0;
	virtual void doEventStep(size_t iEvent) = 0;
	virtual void doPostProcessingSteps() = 0;
	virtual void doReportState(double t) = 0;
};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SYNCHRONOUSAPPROXIMATOR_H_ */
