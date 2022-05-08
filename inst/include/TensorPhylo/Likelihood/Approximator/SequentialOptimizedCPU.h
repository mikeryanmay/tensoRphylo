/*
 * SequentialOptimizedCPU.h
 *
 *  Created on: Nov 19, 2019
 *      Author: xaviermeyer
 *
 *
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SEQUENTIALOPTIMIZEDCPU_H_
#define LIKELIHOOD_APPROXIMATOR_SEQUENTIALOPTIMIZEDCPU_H_

#include "Likelihood/Kernels/CPU/Optimized/EigenIntegrationKernel.hpp"
#include "Likelihood/Kernels/CPU/Optimized/EigenKernels.hpp"
#include "Likelihood/StateTypes/Optimized/EigenState.hpp"
#include "Likelihood/StateTypes/Optimized/EigenStateOperations.hpp"
#include "SynchronousApproximator.h"

#include "Utils/Parallel/Manager.h"

namespace Likelihood {
namespace Approximator {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
class SequentialOptimizedCPU: public SynchronousApproximator {
public:
	SequentialOptimizedCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	   	Conditions::conditionalProbability_t aConditionType,
					Phylogeny::Data::ContainerSharedPtr aPtrData,
					Scheduler::SchedulerSharedPtr aPtrScheduler,
					SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
					Tensor::ContainerSharedPtr aPtrTensorCont);
	~SequentialOptimizedCPU();

	void setDefaultDeltaT(double aDeltaT);
	size_t getTotalNumberOfIntegrationSteps() const;


private:

	const size_t N_MAX_STATE_VECTOR;

	typedef Likelihood::StateType::Optimized::EigenState<conditionalProbType> stateType_t;
	typedef Likelihood::StateType::Optimized::EigenStateOperations<conditionalProbType> operations_t;
	typedef Likelihood::Kernels::CPU::Optimized::EigenKernels<conditionalProbType, withCladoEvents> kernels_t;
	typedef Likelihood::Kernels::CPU::Optimized::EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents> intKernel_t;

	stateType_t probState;
	kernels_t kernels;
	intKernel_t intKernel;

	Likelihood::Integrator::Base<stateType_t, intKernel_t, operations_t>* ptrIntegrator;

	void doPreProcessingSteps();
	void doIntegrationStep(size_t iEdgesLayer);
	void doEventStep(size_t iEvent);
	void doPostProcessingSteps();
	void doReportState(double t);


};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SEQUENTIALOPTIMIZEDCPU_H_ */
