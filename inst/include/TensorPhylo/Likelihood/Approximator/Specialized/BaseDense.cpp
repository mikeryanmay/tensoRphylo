/*
 * BaseDense.cpp
 *
 *  Created on: Apr 16, 2020
 *      Author: meyerx
 */

#include "BaseDense.h"

#include <Eigen/Core>

#include "Utils/Parallel/Manager.h"
#include "Tensor/Container.h"
#include "Utils/MemoryPool/EigenCPU.h"
#include "Likelihood/Scheduler/Event.h"
#include "Likelihood/Scheduler/BaseScheduler.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

BaseDense::BaseDense(processType_t aExtantProcessType,
		Scheduler::SchedulerSharedPtr aPtrScheduler,
		SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		Tensor::ContainerSharedPtr aPtrTensorCont) :
		BaseSpecialized(aExtantProcessType, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont) {

}

BaseDense::~BaseDense() {
}

void BaseDense::defineProbabilities(double time, Eigen::VectorXd &vecProb) {


	assert(ready);

	int iStepper = -1;
	for(size_t iS=0; iS<stepTimes.size(); ++iS) {
		if(time >= stepTimes[iS].first && time <= stepTimes[iS].second) {
			iStepper = iS;
			break;
		}
	}

	assert(iStepper >= 0);
	// Do stuff

	Likelihood::StateType::Vector::EigenState aState;
#ifdef _OPENMP
	if(omp_in_parallel()) {
		#pragma omp critical (steppers)
		{
			steppers[iStepper].calc_state(time, aState,
					stepStates[iStepper].first, stepDerivs[iStepper].first, stepTimes[iStepper].first,
					stepStates[iStepper].second, stepDerivs[iStepper].second, stepTimes[iStepper].second);
		}
	} else {
				steppers[iStepper].calc_state(time, aState,
						stepStates[iStepper].first, stepDerivs[iStepper].first, stepTimes[iStepper].first,
						stepStates[iStepper].second, stepDerivs[iStepper].second, stepTimes[iStepper].second);
	}
#else
		steppers[iStepper].calc_state(time, aState,
				stepStates[iStepper].first, stepDerivs[iStepper].first, stepTimes[iStepper].first,
				stepStates[iStepper].second, stepDerivs[iStepper].second, stepTimes[iStepper].second);
#endif

	vecProb = aState.getStateProb();
}


void BaseDense::doPreprocessingStep() {
	BaseSpecialized::doPreprocessingStep();
	// Clean memory before start
	steppers.clear();
	stepTimes.clear();
	stepStates.clear();
	stepDerivs.clear();
}

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */
