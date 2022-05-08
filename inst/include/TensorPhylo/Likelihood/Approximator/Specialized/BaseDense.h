/*
 * BaseDense.h
 *
 *  Created on: Apr 16, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_BASEDENSE_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_BASEDENSE_H_

#include <Eigen/Core>

#include "BaseSpecialized.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class BaseDense : public BaseSpecialized {
public:
	BaseDense(processType_t aExtantProcessType,
			Scheduler::SchedulerSharedPtr aPtrScheduler,
			SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			Tensor::ContainerSharedPtr aPtrTensorCont);
	virtual ~BaseDense();

	void defineProbabilities(double time, Eigen::VectorXd &vecProb);

protected: // member


	void doPreprocessingStep();
	std::vector<rkd5_stepper_t> steppers;
	std::vector< std::pair<double, double> > stepTimes;
	std::vector< std::pair<StateType, StateType> > stepStates; // X(t), X(t+dt)
	std::vector< std::pair<StateType, StateType> > stepDerivs; // X'(t), X'(t+dt)

};

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_BASEDENSE_H_ */
