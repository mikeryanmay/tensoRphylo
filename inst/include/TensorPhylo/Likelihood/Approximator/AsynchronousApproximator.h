/*
 * AsynchronousApproximator.h
 *
 *  Created on: Apr 17, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_ASYNCHRONOUSAPPROXIMATOR_H_
#define LIKELIHOOD_APPROXIMATOR_ASYNCHRONOUSAPPROXIMATOR_H_

#include "../../Utils/Profiling/CustomProfiling.h"
#include "BaseApproximator.h"

#include "IncFwdLikelihoodApproximator.h"
#include "Likelihood/Scheduler/DAG/IncFwdDAG.h"

namespace Likelihood {
namespace Approximator {

class AsynchronousApproximator: public BaseApproximator {
public:
	AsynchronousApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
			 Conditions::conditionalProbability_t aConditionType,
			 Phylogeny::Data::ContainerSharedPtr aPtrData,
			 Scheduler::SchedulerSharedPtr aPtrScheduler,
			 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			 Tensor::ContainerSharedPtr aPtrTensorCont);
	virtual ~AsynchronousApproximator();


	double approximateLogLikelihood();

protected:

	Specialized::SpecializedApproxManagerSharedPtr ptrSApproxManager;
	Likelihood::Scheduler::DAG::DAGSharedPtr ptrDAG;
	Likelihood::Scheduler::DAG::SchedulerDAGSharedPtr ptrSchedDAG;


	virtual void doPreProcessingSteps() = 0;
	virtual void doIntegrationStep(Likelihood::Scheduler::DAG::NodeDAG* task) = 0;
	virtual void doEventStep(Likelihood::Scheduler::DAG::NodeDAG* task) = 0;
	virtual void doPostProcessingSteps() = 0;
	virtual void doReportState(Likelihood::Scheduler::DAG::NodeDAG* task) = 0;

};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_ASYNCHRONOUSAPPROXIMATOR_H_ */
