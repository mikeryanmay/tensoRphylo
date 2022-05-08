/*
 * AdaptiveSegmentFW.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVESEGMENTFW_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVESEGMENTFW_H_

#include "BaseSpecialized.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Likelihood/Kernels/CPU/Specialized/EigenIntegrationKernelForward.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class AdaptiveSegmentFW {
public:
	AdaptiveSegmentFW(double aStarTime, double aEndTime,
			 const Likelihood::StateType::Vector::EigenState &aInitState,
			 Tensor::ContainerSharedPtr aPtrTensorCont,
			 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrApproxManager);
	~AdaptiveSegmentFW();

	const Likelihood::StateType::Vector::EigenState& getStateProb() const;

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


	double startTime, endTime;
	StateType localState;
	IntegrationKernelType intKernel;

	void compute();

};

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVESEGMENTFW_H_ */
