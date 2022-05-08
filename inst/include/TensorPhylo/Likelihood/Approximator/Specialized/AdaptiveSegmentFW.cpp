/*
 * AdaptiveSegmentFW.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include <cmath>

#include "Data/Structure/Node.h"
#include "Tensor/Container.h"
#include "../BaseApproximator.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "AdaptiveSegmentFW.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.h"


namespace Likelihood {
namespace Approximator {
namespace Specialized {

AdaptiveSegmentFW::AdaptiveSegmentFW(double aStarTime, double aEndTime,
		 const Likelihood::StateType::Vector::EigenState &aInitState,
		 Tensor::ContainerSharedPtr aPtrTensorCont,
		 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrApproxManager) :
						 startTime(aStarTime), endTime(aEndTime),
						 localState(aInitState),
						 intKernel(aPtrApproxManager, aPtrTensorCont) {


	compute();

}

AdaptiveSegmentFW::~AdaptiveSegmentFW() {
}

void AdaptiveSegmentFW::compute() {
	if(startTime == endTime) return; // nothing to be done

	// Create dense integrator
	Likelihood::Integrator::Base<StateType, IntegrationKernelType, OperationType>* ptrIntegrator;
	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<StateType, IntegrationKernelType, OperationType >(BaseApproximator::DEFAULT_ABS_TOLERANCE,
			BaseApproximator::DEFAULT_REL_TOLERANCE,
			-BaseApproximator::DEFAULT_DELTA_T,
			Integrator::RUNGE_KUTTA_DOPRI5);

	ptrIntegrator->integrate(endTime, startTime, localState, intKernel);

	delete ptrIntegrator;
}


const Likelihood::StateType::Vector::EigenState& AdaptiveSegmentFW::getStateProb() const {
	return localState;
}

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */
