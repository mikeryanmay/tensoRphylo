/*
 * EigenInterationKernelForward.h
 *
 *  Created on: May 10, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_SPECIALIZED_EIGENINTEGRATIONKERNELFORWARD_H_
#define LIKELIHOOD_KERNELS_CPU_SPECIALIZED_EIGENINTEGRATIONKERNELFORWARD_H_

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Data/Reader/IncFwdPhyloReader.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

class EigenIntegrationKernelForward {
public:
	EigenIntegrationKernelForward(
			Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager,
			Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~EigenIntegrationKernelForward();

	void operator() ( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t );

private:

	Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr ptrSApproxManager;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

	void doIntegrationStep( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t );
};

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_SPECIALIZED_EIGENINTEGRATIONKERNELFORWARD_H_ */
