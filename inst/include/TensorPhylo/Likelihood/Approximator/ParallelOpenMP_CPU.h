/*
 * ParallelOpenMP.h
 *
 *  Created on: Nov 19, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_PARALLEL_OPENMP_CPU_H_
#define LIKELIHOOD_APPROXIMATOR_PARALLEL_OPENMP_CPU_H_

#include "Utils/Profiling/CustomProfiling.h"
#include "Likelihood/StateTypes/OpenMP/EigenState.hpp"
#include "Likelihood/StateTypes/OpenMP/EigenStateOperations.hpp"
#include "Likelihood/Kernels/CPU/OpenMP/EigenIntegrationKernel.hpp"
#include "Likelihood/Kernels/CPU/OpenMP/EigenKernels.hpp"
#include "SynchronousApproximator.h"

#if defined(_OPENMP)

namespace Likelihood {
namespace Approximator {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
class ParallelOpenMP: public SynchronousApproximator {
public:
	ParallelOpenMP(Likelihood::Integrator::integrationScheme_t aIntScheme,
			   	   Conditions::conditionalProbability_t aConditionType,
				   Phylogeny::Data::ContainerSharedPtr aPtrData,
				   Scheduler::SchedulerSharedPtr aPtrScheduler,
				   SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
				   Tensor::ContainerSharedPtr aPtrTensorCont);
	~ParallelOpenMP();

	void setDefaultDeltaT(double aDeltaT);
	size_t getTotalNumberOfIntegrationSteps() const;

private:

	const size_t N_MAX_STATE_VECTOR;

	typedef Likelihood::StateType::OpenMP::EigenState<conditionalProbType> stateType_t;
	typedef Likelihood::StateType::OpenMP::EigenStateOperations<conditionalProbType> operations_t;
	typedef Likelihood::Kernels::CPU::OpenMP::EigenKernels<conditionalProbType, withCladoEvents> kernels_t;
	typedef Likelihood::Kernels::CPU::OpenMP::EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents> intKernel_t;

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

#endif //defined(_OPENMP)
#endif /* LIKELIHOOD_APPROXIMATOR_PARALLEL_OPENMP_CPU_H_ */
