/*
 * SpecializedApproxManager.cpp
 *
 *  Created on: Apr 16, 2020
 *      Author: meyerx
 */

#include "../Factory.h"
#include "SpecializedApproxManager.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

SpecializedApproxManager::SpecializedApproxManager(
		Conditions::conditionalProbability_t aConditionType,
		   Scheduler::SchedulerSharedPtr aPtrScheduler,
		   SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		   Tensor::ContainerSharedPtr aPtrTensorCont) :
	conditionType(aConditionType) {
	init(aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);
}

SpecializedApproxManager::~SpecializedApproxManager() {
}

void SpecializedApproxManager::compute() {

#ifdef _OPENMP

	if(Utils::Parallel::Manager::getInstance()->useOpenMP()) {
		#pragma omp parallel
		{ // START OMP PARALLEL
			#pragma omp single
			{ // make sure that all task are only executed once
				#pragma omp task
				{
					// Create dense U
					denseU->compute();
				}

				#pragma omp task
				{
					// Create the other as a function of the conditioning probability
					if(conditionType == Conditions::ROOT_SURVIVAL ||
							conditionType == Conditions::ROOT_MRCA ||
							conditionType == Conditions::STEM_SURVIVAL) {
						specializedUHat->compute();
					}
				}

				#pragma omp task
				{
					if(conditionType == Conditions::STEM_TWO_SAMPLES) {
						specializedPairedUS->compute();
					}
				}

				#pragma omp task
				{
					if(conditionType == Conditions::STEM_TWO_EXT_SAMPLES) {
						specializedPairedUSHat->compute();
					}
				}
			}
		} // END OMP PARALLEL
	} else {
		// Create dense U
		denseU->compute();

		// Create the other as a function of the conditioning probability
		if(conditionType == Conditions::ROOT_SURVIVAL ||
				conditionType == Conditions::ROOT_MRCA ||
				conditionType == Conditions::STEM_SURVIVAL) {
			specializedUHat->compute();
		}

		if(conditionType == Conditions::STEM_TWO_SAMPLES) {
			specializedPairedUS->compute();
		}

		if(conditionType == Conditions::STEM_TWO_EXT_SAMPLES) {
			specializedPairedUSHat->compute();
		}
	}
#else
	// Create dense U
	denseU->compute();

	// Create the other as a function of the conditioning probability
	if(conditionType == Conditions::ROOT_SURVIVAL ||
			conditionType == Conditions::ROOT_MRCA ||
			conditionType == Conditions::STEM_SURVIVAL) {
		specializedUHat->compute();
	}

	if(conditionType == Conditions::STEM_TWO_SAMPLES) {
		specializedPairedUS->compute();
	}

	if(conditionType == Conditions::STEM_TWO_EXT_SAMPLES) {
		specializedPairedUSHat->compute();
	}
#endif
}

DenseUnobservedSharedPtr SpecializedApproxManager::getDenseUnobserved(processType_t aExtantProcessType) {
	assert(aExtantProcessType == FULL_PROCESS);
	return denseU;
}

BaseSpecializedSharedPtr SpecializedApproxManager::getSpecializedUnobserved(processType_t aExtantProcessType) {
	if(aExtantProcessType == EXTANT_PROCESS) {
		if(conditionType == Conditions::ROOT_SURVIVAL ||
				conditionType == Conditions::ROOT_MRCA ||
				conditionType == Conditions::STEM_SURVIVAL) {
			return specializedUHat;
		} else {
			assert(false && "specialized u_hat is not to be computed for this conditioning.");
		}
		return BaseSpecializedSharedPtr();
	} else {
		return denseU;
	}
}

AdaptivePairedUSSharedPtr SpecializedApproxManager::getSpecializedPairedUS(processType_t aExtantProcessType) {
	if(aExtantProcessType == EXTANT_PROCESS) {
		assert(conditionType == Conditions::STEM_TWO_EXT_SAMPLES);
		return specializedPairedUSHat;
	} else {
		assert(conditionType == Conditions::STEM_TWO_SAMPLES);
		return specializedPairedUS;
	}
}

void SpecializedApproxManager::init(Scheduler::SchedulerSharedPtr aPtrScheduler,
							   SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
							   Tensor::ContainerSharedPtr aPtrTensorCont) {

	// Create dense U
	denseU = Factory::createDenseUnobservedCPU(FULL_PROCESS, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);

	// Create the other as a function of the conditioning probability
	if(conditionType == Conditions::ROOT_SURVIVAL ||
			conditionType == Conditions::ROOT_MRCA ||
			conditionType == Conditions::STEM_SURVIVAL) {
		specializedUHat = Factory::createAdaptiveUnobservedCPU(EXTANT_PROCESS, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);
	}

	if(conditionType == Conditions::STEM_TWO_SAMPLES) {
		specializedPairedUS = Factory::createAdaptivePairedUS(FULL_PROCESS, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);
	}

	if(conditionType == Conditions::STEM_TWO_EXT_SAMPLES) {
		specializedPairedUSHat = Factory::createAdaptivePairedUS(EXTANT_PROCESS, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont);
	}


}

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */
