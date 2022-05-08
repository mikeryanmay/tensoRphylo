/*
 * SpecializedApproxManager.h
 *
 *  Created on: Apr 16, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_SPECIALIZEDAPPROXMANAGER_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_SPECIALIZEDAPPROXMANAGER_H_

#include "Tensor/IncFwdTensor.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Likelihood/ConditionTypes/ConditionType.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class SpecializedApproxManager {
public:
	SpecializedApproxManager(Conditions::conditionalProbability_t aConditionType,
					   Scheduler::SchedulerSharedPtr aPtrScheduler,
					   SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
					   Tensor::ContainerSharedPtr aPtrTensorCont);

	~SpecializedApproxManager();

	void compute();

	DenseUnobservedSharedPtr getDenseUnobserved(processType_t aExtantProcessType);

	BaseSpecializedSharedPtr getSpecializedUnobserved(processType_t aExtantProcessType);
	AdaptivePairedUSSharedPtr getSpecializedPairedUS(processType_t aExtantProcessType);

private:

	Conditions::conditionalProbability_t conditionType;

	DenseUnobservedSharedPtr denseU;
	BaseSpecializedSharedPtr specializedUHat;
	AdaptivePairedUSSharedPtr specializedPairedUS, specializedPairedUSHat;

	void init(Scheduler::SchedulerSharedPtr aPtrScheduler,
			   SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			   Tensor::ContainerSharedPtr aPtrTensorCont);


};

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_SPECIALIZEDAPPROXMANAGER_H_ */
