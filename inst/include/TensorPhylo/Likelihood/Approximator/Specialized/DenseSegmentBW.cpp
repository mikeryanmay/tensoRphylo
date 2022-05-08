/*
 * DenseSegmentBW.cpp
 *
 *  Created on: May 7, 2020
 *      Author: meyerx
 */

#include "DenseSegmentBW.h"

#include "../BaseApproximator.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Utils/Parallel/Manager.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

DenseSegmentBW::DenseSegmentBW(
		 double aStarTime, double aEndTime,
		 const Likelihood::StateType::Vector::EigenState &aInitState,
		 Tensor::ContainerSharedPtr aPtrTensorCont,
		 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager) :
				 startTime(aStarTime), endTime(aEndTime),
				 localState(aInitState),
				 intKernel(aPtrSApproxManager, aPtrTensorCont) {
	compute();

}

DenseSegmentBW::~DenseSegmentBW() {
}

void DenseSegmentBW::compute() {
	if(startTime == endTime) return; // nothing to be done

	// Create dense integrator
	Likelihood::Integrator::Base<StateType, IntegrationKernelType, OperationType>* ptrIntegrator;
	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<StateType, IntegrationKernelType, OperationType >(BaseApproximator::DEFAULT_ABS_TOLERANCE,
			BaseApproximator::DEFAULT_REL_TOLERANCE,
			BaseApproximator::DEFAULT_DELTA_T,
			Integrator::DENSE_RUNGE_KUTTA_DOPRI5);

	// Create dense integrator
	Likelihood::Integrator::DenseRungeKuttaDOPRI5<StateType, IntegrationKernelType, OperationType> *ptrDRK =
			dynamic_cast< Likelihood::Integrator::DenseRungeKuttaDOPRI5<StateType, IntegrationKernelType, OperationType>* >(ptrIntegrator);
	assert(ptrDRK != NULL);

	// Copy initial state and integrate
	ptrDRK->integrate(startTime, endTime, localState, intKernel);
	// Check that we end up to the same state prob as in the forward pass

	ptrDRK->transferStepsMemory(steppers, stepTimes, stepStates, stepDerivs);

	delete ptrIntegrator;
}


void DenseSegmentBW::defineProbabilities(double time, Eigen::VectorXd &vecProb) {

	assert(time >= stepTimes.front().first && time <= stepTimes.back().second);

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

const Likelihood::StateType::Vector::EigenState& DenseSegmentBW::getStateProb() const {
	return localState;
}


} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */
