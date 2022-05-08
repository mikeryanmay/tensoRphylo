/*
 * SingleStepSegmentFW.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include <cmath>

#include "Data/Structure/Node.h"
#include "Tensor/Container.h"
#include "../BaseApproximator.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "SingleStepSegmentFW.h"

#include <boost/numeric/odeint/stepper/generation/make_dense_output.hpp>
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.h"
#include "Likelihood/CustomIntegrators/AdaptiveIntegrators.hpp"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

SingleStepSegmentFW::SingleStepSegmentFW(double aDT,
		 Tensor::ContainerSharedPtr aPtrTensorCont,
		 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrApproxManager) :
						 intKernel(aPtrApproxManager, aPtrTensorCont),
						 denseSteppper(boost::numeric::odeint::make_dense_output( 1E-7 , 1E-7 , rkd5_stepper_t())) {

	first = true;
	curTStart = curTEnd = 0.;
	dt = aDT;

}

SingleStepSegmentFW::~SingleStepSegmentFW() {
}

double SingleStepSegmentFW::doNextStep(Likelihood::StateType::Vector::EigenState& intOutState, double &startTime, const double &endTime) {

	curTStart = startTime;

	// Do one dt step
	curTEnd = integrate_adaptive_stepwise( boost::ref(denseSteppper), boost::ref(intKernel), intOutState, startTime, endTime , dt, first);

	first = false;

	return curTEnd;
}

void SingleStepSegmentFW::reset() {
	first = true;
	curTStart = curTEnd = 0.;
}


void SingleStepSegmentFW::defineProbabilities(double time, Likelihood::StateType::Vector::EigenState& outState) {

	assert(curTStart >= time && curTEnd <= time && !first);
	denseSteppper.calc_state(time, outState);

}


} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */
