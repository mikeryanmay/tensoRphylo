/*
 * AdaptiveUnobserved.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVEOUTPUTUNOBSERVED_H_
#define LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVEOUTPUTUNOBSERVED_H_

#include "BaseSpecialized.h"
#include "Likelihood/Kernels/CPU/Specialized/EigenKernelsUnobserved.h"
#include "Likelihood/Kernels/CPU/Specialized/IntegrationKernelUnobserved.h"

namespace Likelihood {
namespace Approximator {
namespace Specialized {

class AdaptiveUnobserved : public BaseSpecialized {
public:
	AdaptiveUnobserved(processType_t aExtantProcessType,
					Scheduler::SchedulerSharedPtr aPtrScheduler,
					SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
					Tensor::ContainerSharedPtr aPtrTensorCont);
	~AdaptiveUnobserved();

private:

	// Kernels
	typedef Likelihood::Kernels::CPU::Specialized::EigenKernelsUnobserved KernelsType;
	typedef Likelihood::Kernels::CPU::Specialized::IntegrationKernelUnobserved IntegrationKernelsType;

	KernelsType kernels;
	IntegrationKernelsType intKernel;

	void doIntegrationStep(double startTime, double endTime);
	void doEventStep(size_t iEvent);

};

} /* namespace Specialized */
} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SPECIALIZED_ADAPTIVEOUTPUTUNOBSERVED_H_ */
