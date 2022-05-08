/*
 * EigenKernelsSingleton.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_EIGENKERNELSSINGLETON_H_
#define LIKELIHOOD_KERNELS_CPU_EIGENKERNELSSINGLETON_H_

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

class EigenKernelsSingleton {
public:
	EigenKernelsSingleton(Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
					Likelihood::Approximator::Specialized::DenseUnobservedSharedPtr aPtrUnobsProb,
					SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
					Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~EigenKernelsSingleton();

	void setInitialCondition(Likelihood::StateType::Vector::EigenState &x);
	void computeMassSpeciation(double t, Likelihood::StateType::Vector::EigenState &x);
	void computeMassExtinction(double t, Likelihood::StateType::Vector::EigenState &x);
	void computeMassSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x);
	void computeMassDestrSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x);

private:

	const Likelihood::Approximator::Specialized::processType_t PROCESS_TYPE;

	Likelihood::Approximator::Specialized::DenseUnobservedSharedPtr ptrDenseUnobs;
	SynchronousEvents::ContainerSharedPtr ptrSynchEventContainer;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

};

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_EIGENKERNELSSINGLETON_H_ */
