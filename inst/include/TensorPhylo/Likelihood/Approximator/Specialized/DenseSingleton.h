/*
 * AdaptiveSingleton.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_DENSEOUTPUTSINGLETON_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_DENSEOUTPUTSINGLETON_H_

#include "BaseDense.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "Likelihood/Kernels/CPU/Specialized/EigenKernelsSingleton.h"
#include "Likelihood/Kernels/CPU/Specialized/IntegrationKernelSingleton.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class DenseSingleton : public BaseDense {
public:
	DenseSingleton(processType_t aExtantProcessType,
					Scheduler::SchedulerSharedPtr aPtrScheduler,
					Likelihood::Approximator::Specialized::DenseUnobservedSharedPtr aPtrDenseUnobs,
					SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
					Tensor::ContainerSharedPtr aPtrTensorCont);
	~DenseSingleton();

private:

	DenseUnobservedSharedPtr ptrDenseUnobs;

	// State
	typedef Likelihood::Kernels::CPU::Specialized::EigenKernelsSingleton KernelsType;
	typedef Likelihood::Kernels::CPU::Specialized::IntegrationKernelSingleton IntegrationKernelType;

	KernelsType kernels;
	IntegrationKernelType intKernel;

	void doIntegrationStep(double startTime, double endTime);
	void doEventStep(size_t iEvent);

};

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_DENSEOUTPUTSINGLETON_H_ */
