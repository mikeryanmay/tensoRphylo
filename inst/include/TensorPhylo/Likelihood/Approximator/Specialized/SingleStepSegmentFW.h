/*
 * SingleStepSegmentFW.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_SINGLESTEPSEGMENTFW_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_SINGLESTEPSEGMENTFW_H_

#include "BaseSpecialized.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Likelihood/Kernels/CPU/Specialized/EigenIntegrationKernelForward.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class SingleStepSegmentFW {
public:
	SingleStepSegmentFW(double aDT,
			 Tensor::ContainerSharedPtr aPtrTensorCont,
			 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrApproxManager);
	~SingleStepSegmentFW();

	double doNextStep(Likelihood::StateType::Vector::EigenState& intOutState, double &startTime, const double &endTime);
	void reset();

	void defineProbabilities(double time, Likelihood::StateType::Vector::EigenState& outState);

protected: // member

	// State
	typedef Likelihood::StateType::Vector::EigenState StateType;
	typedef Likelihood::StateType::Vector::EigenStateOperations OperationType;
	typedef Likelihood::Kernels::CPU::Specialized::EigenIntegrationKernelForward IntegrationKernelType;

	// Stepper
	typedef boost::numeric::odeint::runge_kutta_dopri5<
			StateType ,
			double ,
			StateType,
			double ,
			boost::numeric::odeint::vector_space_algebra,
			OperationType,
			boost::numeric::odeint::always_resizer > rkd5_stepper_t;

	typedef boost::numeric::odeint::controlled_runge_kutta< rkd5_stepper_t > controlledRDK5_t;


	typedef boost::numeric::odeint::dense_output_runge_kutta< controlledRDK5_t > denseRDK5_t;

	bool first;
	double curTStart, curTEnd, dt;

	IntegrationKernelType intKernel;
	denseRDK5_t denseSteppper;

};

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_SINGLESTEPSEGMENTFW_H_ */
