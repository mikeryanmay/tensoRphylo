/*
 * BranchwiseCPU.h
 *
 *  Created on: April 17, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SEQUENTIALBRANCHWISECPU_H_
#define LIKELIHOOD_APPROXIMATOR_SEQUENTIALBRANCHWISECPU_H_

#include "AsynchronousApproximator.h"
#include "IncFwdLikelihoodApproximator.h"
#include "Likelihood/Scheduler/DAG/IncFwdDAG.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/Kernels/CPU/IncEigenKernels.h"
#include "Likelihood/StateTypes/Vector/EigenStateOperations.hpp"

namespace Likelihood {
namespace Approximator {

class SequentialBranchwiseCPU: public AsynchronousApproximator {
public:
	SequentialBranchwiseCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	  Conditions::conditionalProbability_t aConditionType,
				  Phylogeny::Data::ContainerSharedPtr aPtrData,
				  Scheduler::SchedulerSharedPtr aPtrScheduler,
				  SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
				  Tensor::ContainerSharedPtr aPtrTensorCont);
	~SequentialBranchwiseCPU();

	void setDefaultDeltaT(double aDeltaT);
	size_t getTotalNumberOfIntegrationSteps() const;


private:

	typedef Likelihood::StateType::Vector::EigenState stateType_t;
	typedef Likelihood::Kernels::CPU::Branchwise::EigenKernels kernels_t;
	typedef Likelihood::Kernels::CPU::Branchwise::EigenIntegrationKernel intKernel_t;
	typedef Likelihood::StateType::Vector::EigenStateOperations operations_t;

	stateType_t probState;

	kernels_t kernels;
	intKernel_t intKernel;

	Likelihood::Integrator::Base<stateType_t, intKernel_t, operations_t>* ptrIntegrator;

	void doPreProcessingSteps();
	void doIntegrationStep(Likelihood::Scheduler::DAG::NodeDAG* task);
	void doEventStep(Likelihood::Scheduler::DAG::NodeDAG* task) ;
	void doPostProcessingSteps();
	void doReportState(Likelihood::Scheduler::DAG::NodeDAG* task);

};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SEQUENTIALBRANCHWISECPU_H_ */
