/*
 * ParallelAsynchApproximator.h
 *
 *  Created on: Apr 17, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_PARALLELASYNCHAPPROXIMATOR_H_
#define LIKELIHOOD_APPROXIMATOR_PARALLELASYNCHAPPROXIMATOR_H_

#include "BaseApproximator.h"

#include "Utils/Parallel/Manager.h"
#if defined(_OPENMP)

#include "IncFwdLikelihoodApproximator.h"
#include "Likelihood/Scheduler/DAG/IncFwdDAG.h"

namespace Likelihood {
namespace Approximator {

class ParallelAsynchApproximator: public BaseApproximator {
public:
	ParallelAsynchApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
			 Conditions::conditionalProbability_t aConditionType,
			 Phylogeny::Data::ContainerSharedPtr aPtrData,
			 Scheduler::SchedulerSharedPtr aPtrScheduler,
			 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			 Tensor::ContainerSharedPtr aPtrTensorCont);
	virtual ~ParallelAsynchApproximator();


	double approximateLogLikelihood();

protected:

	Specialized::SpecializedApproxManagerSharedPtr ptrSApproxManager;
	Likelihood::Scheduler::DAG::DAGSharedPtr ptrDAG;
	Likelihood::Scheduler::DAG::ParallelSchedulerDAGSharedPtr ptrSchedDAG;


	virtual void doPreProcessingSteps() = 0;
	virtual void doIntegrationStep(Likelihood::Scheduler::DAG::NodeDAG* task) = 0;
	virtual void doEventStep(Likelihood::Scheduler::DAG::NodeDAG* task) = 0;
	virtual void doPostProcessingSteps() = 0;
	virtual void doReportState(Likelihood::Scheduler::DAG::NodeDAG* task) = 0;

};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif

#endif /* LIKELIHOOD_APPROXIMATOR_PARALLELASYNCHAPPROXIMATOR_H_ */
