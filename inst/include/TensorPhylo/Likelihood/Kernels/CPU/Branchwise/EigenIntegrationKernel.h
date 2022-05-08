/*
 * EigenInterationKernel.h
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_BRANCHWISE_EIGENINTEGRATIONKERNEL_H_
#define LIKELIHOOD_KERNELS_CPU_BRANCHWISE_EIGENINTEGRATIONKERNEL_H_

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
namespace Branchwise {

class EigenIntegrationKernel {
public:
	EigenIntegrationKernel(Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager,
						   Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~EigenIntegrationKernel();

	void operator() ( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t );

	std::pair<double, bool> getStiffnessRatio(double t, const Eigen::VectorXd &u) const;

private:

	Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr ptrSApproxManager;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

	void doIntegrationStep( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t );
};

} /* namespace Branchwise */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_BRANCHWISE_EIGENINTEGRATIONKERNEL_H_ */
