/*
 * EigenKernelsUnobserved.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_EIGENKERNELSUNOBSERVED_H_
#define LIKELIHOOD_KERNELS_CPU_EIGENKERNELSUNOBSERVED_H_

#include "Tensor/IncFwdTensor.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"


namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

class EigenKernelsUnobserved {
public:
	EigenKernelsUnobserved(
					Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
					SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
					Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~EigenKernelsUnobserved();

	void setInitialCondition(Likelihood::StateType::Vector::EigenState &x);
	void computeMassSpeciation(double t, Likelihood::StateType::Vector::EigenState &x);
	void computeMassExtinction(double t, Likelihood::StateType::Vector::EigenState &x);
	void computeMassSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x);
	void computeMassDestrSamplingEvent(double t, Likelihood::StateType::Vector::EigenState &x);

private:

	const Likelihood::Approximator::Specialized::processType_t PROCESS_TYPE;

	SynchronousEvents::ContainerSharedPtr ptrSynchEventContainer;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

};

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_EIGENKERNELSUNOBSERVED_H_ */
