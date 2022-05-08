/*
 * Factory.cpp
 *
 *  Created on: Dec 2, 2019
 *      Author: meyerx
 */

#include "Factory.h"

#include "AutoTuningApproximator.h"
#include "ParallelBranchwiseCPU.h"
#include "SequentialBranchwiseCPU.h"
#include "SequentialOptimizedCPU.hpp"
#include "Specialized/AdaptiveUnobserved.h"
#include "Specialized/AdaptiveSingleton.h"
#include "Specialized/DenseUnobserved.h"
#include "Specialized/DenseSingleton.h"
#include "StochasticMapping/SequentialStochasticMappingCPU.h"

#if defined(_OPENMP)
#include "ParallelOpenMP_CPU.hpp"
#endif

namespace Likelihood {
namespace Approximator {

Factory::Factory() {

}

Factory::~Factory() {
}

ApproximatorSharedPtr Factory::createSequentialTemplateCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
													 Conditions::conditionalProbability_t aConditionType,
													 Phylogeny::Data::ContainerSharedPtr aPtrData,
													 Scheduler::SchedulerSharedPtr aPtrScheduler,
													 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
													 Tensor::ContainerSharedPtr aPtrTensorCont) {
	return optimizedApproximatorFactory(aIntScheme,  aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);
}

ApproximatorSharedPtr Factory::createSequentialBranchwiseCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
	   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
													 Phylogeny::Data::ContainerSharedPtr aPtrData,
													 Scheduler::SchedulerSharedPtr aPtrScheduler,
													 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
													 Tensor::ContainerSharedPtr aPtrTensorCont) {
	ApproximatorSharedPtr ptr(new Likelihood::Approximator::SequentialBranchwiseCPU(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}


#if defined(_OPENMP)
ApproximatorSharedPtr Factory::createParallelOpenMPCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
												     Conditions::conditionalProbability_t aConditionType,
													 Phylogeny::Data::ContainerSharedPtr aPtrData,
													 Scheduler::SchedulerSharedPtr aPtrScheduler,
													 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
													 Tensor::ContainerSharedPtr aPtrTensorCont){
	return openmpApproximatorFactory(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);
}

ApproximatorSharedPtr Factory::createParallelBranchwiseCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
	   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
													 Phylogeny::Data::ContainerSharedPtr aPtrData,
													 Scheduler::SchedulerSharedPtr aPtrScheduler,
													 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
													 Tensor::ContainerSharedPtr aPtrTensorCont) {
	ApproximatorSharedPtr ptr(new Likelihood::Approximator::ParallelBranchwiseCPU(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

#endif

ApproximatorSharedPtr Factory::createAutoTuningCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
	   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
													 Phylogeny::Data::ContainerSharedPtr aPtrData,
													 Scheduler::SchedulerSharedPtr aPtrScheduler,
													 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
													 Tensor::ContainerSharedPtr aPtrTensorCont) {
	ApproximatorSharedPtr ptr(new Likelihood::Approximator::AutoTuningApproximator(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

DenseUnobsApproximatorSharedPtr Factory::createDenseUnobservedCPU(Specialized::processType_t aExtantProcessType,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont) {
	DenseUnobsApproximatorSharedPtr ptr(new Likelihood::Approximator::Specialized::DenseUnobserved(aExtantProcessType, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

DenseApproximatorSharedPtr Factory::createDenseSingletonCPU(Specialized::processType_t aExtantProcessType,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 DenseUnobsApproximatorSharedPtr aPtrDenseUnobsApprox,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont) {
	DenseApproximatorSharedPtr ptr(new Likelihood::Approximator::Specialized::DenseSingleton(aExtantProcessType, aPtrScheduler, aPtrDenseUnobsApprox, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

SpecializedApproximatorSharedPtr Factory::createAdaptiveUnobservedCPU(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont) {
	SpecializedApproximatorSharedPtr ptr(new Likelihood::Approximator::Specialized::AdaptiveUnobserved(aExtantProcessType, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

SpecializedApproximatorSharedPtr Factory::createAdaptiveSingletonCPU(Specialized::processType_t aExtantProcessType,
															 Scheduler::SchedulerSharedPtr aPtrScheduler,
															 DenseUnobsApproximatorSharedPtr aPtrDenseUnobsApprox,
															 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
															 Tensor::ContainerSharedPtr aPtrTensorCont) {
	SpecializedApproximatorSharedPtr ptr(new Likelihood::Approximator::Specialized::AdaptiveSingleton(aExtantProcessType, aPtrScheduler, aPtrDenseUnobsApprox, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

AdaptivePairedUSSharedPtr Factory::createAdaptivePairedUS(Specialized::processType_t aExtantProcessType,
														 Scheduler::SchedulerSharedPtr aPtrScheduler,
														 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
														 Tensor::ContainerSharedPtr aPtrTensorCont) {
	AdaptivePairedUSSharedPtr ptr(new Likelihood::Approximator::Specialized::AdaptivePairedUS(aExtantProcessType, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}

StochasticMapping::StochasticMappingApproxSharedPtr Factory::createStochasticMappingApprox(Likelihood::Integrator::integrationScheme_t aIntScheme,
	   	   	   	   	   	   	   	   	   	   	   	   	 Conditions::conditionalProbability_t aConditionType,
													 Phylogeny::Data::ContainerSharedPtr aPtrData,
													 Scheduler::SchedulerSharedPtr aPtrScheduler,
													 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
													 Tensor::ContainerSharedPtr aPtrTensorCont) {
	StochasticMapping::StochasticMappingApproxSharedPtr ptr(
			new Likelihood::Approximator::StochasticMapping::SequentialStochasticMappingCPU(
					aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont));
	return ptr;
}


ApproximatorSharedPtr Factory::optimizedApproximatorFactory(
		Likelihood::Integrator::integrationScheme_t aIntScheme,
		Likelihood::Conditions::conditionalProbability_t condProbType,
		Phylogeny::Data::ContainerSharedPtr aPtrData,
		Likelihood::Scheduler::SchedulerSharedPtr scheduler,
		SynchronousEvents::ContainerSharedPtr syncEvents,
		Tensor::ContainerSharedPtr tensors) {

	ApproximatorSharedPtr approx;

	const Tensor::etaStructure_t etaStructureType = tensors->getEtaStructureType();
	const bool withCladoEvents = tensors->getOmega()->getDimensions()[0] > 0;

	if(etaStructureType == Tensor::ETA_SPARSE && withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset(new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::TIME, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_SPARSE && !withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::TIME, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_DENSE && withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::TIME, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_DENSE && !withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::TIME, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_QUASSE && withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::TIME, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_QUASSE && !withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::TIME, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::SequentialOptimizedCPU< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	}

	assert(approx != NULL);
	return approx;
}


#if defined(_OPENMP)
ApproximatorSharedPtr Factory::openmpApproximatorFactory(
		Likelihood::Integrator::integrationScheme_t aIntScheme,
		Likelihood::Conditions::conditionalProbability_t condProbType,
		Phylogeny::Data::ContainerSharedPtr aPtrData,
		Likelihood::Scheduler::SchedulerSharedPtr scheduler,
		SynchronousEvents::ContainerSharedPtr syncEvents,
		Tensor::ContainerSharedPtr tensors) {

	ApproximatorSharedPtr approx;

	const Tensor::etaStructure_t etaStructureType = tensors->getEtaStructureType();
	const bool withCladoEvents = tensors->getOmega()->getDimensions()[0] > 0;

	if(etaStructureType == Tensor::ETA_SPARSE && withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset(new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::TIME, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_SPARSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_SPARSE && !withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::TIME, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_SPARSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_DENSE && withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::TIME, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_DENSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_DENSE && !withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::TIME, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_DENSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_QUASSE && withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::TIME, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_QUASSE, true >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	} else if(etaStructureType == Tensor::ETA_QUASSE && !withCladoEvents) {
		if(condProbType == Likelihood::Conditions::TIME) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::TIME, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SURVIVAL, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_MRCA, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_SURVIVAL) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_SURVIVAL, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_ONE_SAMPLE, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_EXT_SAMPLES, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::STEM_TWO_SAMPLES, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		} else if(condProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA) {
			approx.reset( new Likelihood::Approximator::ParallelOpenMP< Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA, Tensor::ETA_QUASSE, false >(aIntScheme, condProbType, aPtrData, scheduler, syncEvents, tensors));
		}
	}

	assert(approx != NULL);
	return approx;
}
#endif

} /* namespace Approximator */
} /* namespace Likelihood */
