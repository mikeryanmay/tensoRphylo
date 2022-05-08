/*
 * BaseSpecialized.h
 *
 *  Created on: Apr 16, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_BASESPECIALIZED_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_BASESPECIALIZED_H_

#include <Eigen/Core>

#include "Tensor/IncFwdTensor.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/StateTypes/Vector/EigenStateOperations.hpp"

namespace Likelihood {
namespace Scheduler {
namespace DAG {
class NodeDAG;
}
}
}

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class BaseSpecialized {
public:
	BaseSpecialized(
			processType_t aExtantProcessType,
			Scheduler::SchedulerSharedPtr aPtrScheduler,
			SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			Tensor::ContainerSharedPtr aPtrTensorCont);
	virtual ~BaseSpecialized();

	bool isExtantProcess() const;
	bool isFullProcess() const;

	processType_t getProcessType() const;

	void compute();

	const Eigen::VectorXd& getFinalProbability() const;

protected: // type

	// State
	typedef Likelihood::StateType::Vector::EigenState StateType;
	typedef Likelihood::StateType::Vector::EigenStateOperations OperationType;

	// Stepper
	typedef boost::numeric::odeint::runge_kutta_dopri5<
			StateType ,
			double ,
			StateType,
			double ,
			boost::numeric::odeint::vector_space_algebra,
			OperationType,
			boost::numeric::odeint::always_resizer > rkd5_stepper_t;

protected: // member

	const processType_t PROCESS_TYPE;

	bool ready;
	Scheduler::SchedulerSharedPtr ptrScheduler;
	StateType probState;

	virtual void doPreprocessingStep();
	virtual void doIntegrationStep(double startTime, double endTime) = 0;
	virtual void doEventStep(size_t iEvent) = 0;


};

} /* namespace Specialized */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_BASESPECIALIZED_H_ */
