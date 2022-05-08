/*
 * AdaptivePairedUS.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVEPAIREDUS_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVEPAIREDUS_H_

#include "BaseSpecialized.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Likelihood/StateTypes/Matrix/EigenState.h"
#include "Likelihood/StateTypes/Matrix/EigenStateOperations.hpp"
#include "Likelihood/Kernels/CPU/Specialized/EigenKernelsPairedUS.h"
#include "Likelihood/Kernels/CPU/Specialized/IntegrationKernelPairedUS.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class AdaptivePairedUS : public BaseSpecialized {
public:
	AdaptivePairedUS(processType_t aExtantProcessType,
					Scheduler::SchedulerSharedPtr aPtrScheduler,
					SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
					Tensor::ContainerSharedPtr aPtrTensorCont);
	~AdaptivePairedUS();

	const Eigen::VectorXd& getU() const;
	const Eigen::VectorXd& getS() const;

private:

	// Shadowing base on purpose
	// State
	typedef Likelihood::StateType::Matrix::EigenState StateType;
	typedef Likelihood::StateType::Matrix::EigenStateOperations OperationType;
	StateType matrixProbState;

	// Stepper
	typedef boost::numeric::odeint::runge_kutta_dopri5<
			StateType ,
			double ,
			StateType,
			double ,
			boost::numeric::odeint::vector_space_algebra,
			OperationType,
			boost::numeric::odeint::always_resizer > rkd5_stepper_t;

	// State
	typedef Likelihood::Kernels::CPU::Specialized::EigenKernelsPairedUS KernelsType;
	typedef Likelihood::Kernels::CPU::Specialized::IntegrationKernelPairedUS IntegrationKernelType;

	KernelsType kernels;
	IntegrationKernelType intKernel;


	void doPreprocessingStep();
	void doIntegrationStep(double startTime, double endTime);
	void doEventStep(size_t iEvent);

};

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVEPAIREDUS_H_ */
