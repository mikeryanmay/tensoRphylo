/*
 * DenseSegmentBW.h
 *
 *  Created on: May 7, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_DenseSegmentBWBW_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_DenseSegmentBWBW_H_

#include <vector>

#include "Tensor/IncFwdTensor.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/StateTypes/Vector/EigenStateOperations.hpp"
#include "Likelihood/Kernels/CPU/Branchwise/EigenIntegrationKernel.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class DenseSegmentBW {
public:
	DenseSegmentBW(double aStarTime, double aEndTime,
				 const Likelihood::StateType::Vector::EigenState &aInitState,
				 Tensor::ContainerSharedPtr aPtrTensorCont,
				 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager);
	~DenseSegmentBW();

	void defineProbabilities(double time, Eigen::VectorXd &vecProb);
	const Likelihood::StateType::Vector::EigenState& getStateProb() const;

protected: // member

	// State
	typedef Likelihood::StateType::Vector::EigenState StateType;
	typedef Likelihood::StateType::Vector::EigenStateOperations OperationType;
	typedef Likelihood::Kernels::CPU::Branchwise::EigenIntegrationKernel IntegrationKernelType;

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

	std::vector<rkd5_stepper_t> steppers;
	std::vector< std::pair<double, double> > stepTimes;
	std::vector< std::pair<StateType, StateType> > stepStates; // X(t), X(t+dt)
	std::vector< std::pair<StateType, StateType> > stepDerivs; // X'(t), X'(t+dt)

	void compute();

};

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_DenseSegmentBWBW_H_ */
