/*
 * EigenInterationKernel.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_INTEGRATIONKERNELUNOBSERVED_H_
#define LIKELIHOOD_KERNELS_CPU_INTEGRATIONKERNELUNOBSERVED_H_

#include "Tensor/IncFwdTensor.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

class IntegrationKernelUnobserved {
public:
	IntegrationKernelUnobserved(Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
								Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~IntegrationKernelUnobserved();

	void operator() ( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t );

private:

	const Likelihood::Approximator::Specialized::processType_t PROCESS_TYPE;

	Eigen::MatrixXd resFirstContractionU;

	Tensor::ContainerSharedPtr ptrTensorsContainer;

	void doIntegrationStep( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t );
};

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_INTEGRATIONKERNELUNOBSERVED_H_ */
