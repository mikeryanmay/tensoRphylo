/*
 * Factory.h
 *
 *  Created on: Dec 2, 2019
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_FACTORY_H_
#define LIKELIHOOD_APPROXIMATOR_FACTORY_H_

#include "BaseApproximator.h"
#include "Specialized/AdaptivePairedUS.h"
#include "Specialized/BaseSpecialized.h"
#include "Specialized/BaseDense.h"
#include "Specialized/DenseUnobserved.h"

#include "Utils/Parallel/Manager.h"

namespace Likelihood {
namespace Approximator {

typedef boost::shared_ptr<Specialized::BaseSpecialized> SpecializedApproximatorSharedPtr;
typedef boost::shared_ptr<Specialized::BaseDense> DenseApproximatorSharedPtr;
typedef boost::shared_ptr<Specialized::DenseUnobserved> DenseUnobsApproximatorSharedPtr;
typedef boost::shared_ptr<Specialized::AdaptivePairedUS> AdaptivePairedUSSharedPtr;
typedef boost::shared_ptr<BaseApproximator> ApproximatorSharedPtr;

class Factory {
public:

	static ApproximatorSharedPtr createSequentialTemplateCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
														 Phylogeny::Data::ContainerSharedPtr aPtrData,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont);

	static ApproximatorSharedPtr createSequentialBranchwiseCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
														 Phylogeny::Data::ContainerSharedPtr aPtrData,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont);

#if defined(_OPENMP)
	static ApproximatorSharedPtr createParallelOpenMPCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
														 Phylogeny::Data::ContainerSharedPtr aPtrData,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont);


	static ApproximatorSharedPtr createParallelBranchwiseCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
														 Phylogeny::Data::ContainerSharedPtr aPtrData,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont);

#endif

	static ApproximatorSharedPtr createAutoTuningCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
														 Phylogeny::Data::ContainerSharedPtr aPtrData,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont);

	static DenseUnobsApproximatorSharedPtr createDenseUnobservedCPU(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont);

	static DenseApproximatorSharedPtr createDenseSingletonCPU(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 DenseUnobsApproximatorSharedPtr aPtrDenseUnobsApprox,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont);

	static SpecializedApproximatorSharedPtr createAdaptiveUnobservedCPU(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont);

	static SpecializedApproximatorSharedPtr createAdaptiveSingletonCPU(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 DenseUnobsApproximatorSharedPtr aPtrDenseUnobsApprox,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont);

	static AdaptivePairedUSSharedPtr createAdaptivePairedUS(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont);

	static StochasticMapping::StochasticMappingApproxSharedPtr createStochasticMappingApprox(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	   	   	   	   	   	   	   	   	   	 	 Conditions::conditionalProbability_t aConditionType,
															 Phylogeny::Data::ContainerSharedPtr aPtrData,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont);

private:

	Factory();
	~Factory();

	static ApproximatorSharedPtr optimizedApproximatorFactory(
			Likelihood::Integrator::integrationScheme_t aIntScheme,
			Likelihood::Conditions::conditionalProbability_t condProbType,
			Phylogeny::Data::ContainerSharedPtr aPtrData,
			Likelihood::Scheduler::SchedulerSharedPtr scheduler,
			SynchronousEvents::ContainerSharedPtr syncEvents,
			Tensor::ContainerSharedPtr tensors);

#if defined(_OPENMP)
	static ApproximatorSharedPtr openmpApproximatorFactory(
			Likelihood::Integrator::integrationScheme_t aIntScheme,
			Likelihood::Conditions::conditionalProbability_t condProbType,
			Phylogeny::Data::ContainerSharedPtr aPtrData,
			Likelihood::Scheduler::SchedulerSharedPtr scheduler,
			SynchronousEvents::ContainerSharedPtr syncEvents,
			Tensor::ContainerSharedPtr tensors);
#endif

};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_FACTORY_H_ */
