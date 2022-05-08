/*
 * SequentialStochasticMappingCPU.h
 *
 *  Created on: may 6, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_SEQUENTIALSTOCHASTICMAPPINGCPU_H_
#define LIKELIHOOD_APPROXIMATOR_SEQUENTIALSTOCHASTICMAPPINGCPU_H_

#include <unordered_map>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "../../../Utils/Profiling/CustomProfiling.h"
#include "IncFwdStochasticMapping.h"
#include "../AsynchronousApproximator.h"
#include "../IncFwdLikelihoodApproximator.h"
#include "BaseStochasticMappingApproximator.h"
#include "Likelihood/Scheduler/DAG/IncFwdDAG.h"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/Kernels/CPU/IncEigenKernels.h"
#include "Likelihood/StateTypes/Vector/EigenStateOperations.hpp"
#include "Likelihood/Kernels/CPU/Specialized/EigenIntegrationKernelForward.h"

namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

class SequentialStochasticMappingCPU: public AsynchronousApproximator, public BaseStochasticMappingApproximator {
public:
	SequentialStochasticMappingCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
		   	   	  Conditions::conditionalProbability_t aConditionType,
				  Phylogeny::Data::ContainerSharedPtr aPtrData,
				  Scheduler::SchedulerSharedPtr aPtrScheduler,
				  SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
				  Tensor::ContainerSharedPtr aPtrTensorCont);
	~SequentialStochasticMappingCPU();

	void setAlgorithm(stochasticMappingAlgo_t aAlgo);

	void setDefaultDeltaT(double aDeltaT);
	size_t getTotalNumberOfIntegrationSteps() const;

	mapHistories_t drawHistory();
	mapHistories_t drawHistoryAndComputeRates(std::vector<double>& averageLambda, std::vector<double>& averageMu, std::vector<double>& averagePhi, std::vector<double>& averageDelta, std::vector<long>& numChanges);
	mapHistories_t drawAncestralStates();

private:

	typedef Likelihood::StateType::Vector::EigenState stateType_t;
	typedef Likelihood::Kernels::CPU::Branchwise::EigenKernels kernelsBW_t;
	typedef Likelihood::Kernels::CPU::Branchwise::EigenIntegrationKernel intKernelBW_t;
	typedef Likelihood::Kernels::CPU::Specialized::EigenIntegrationKernelForward intKernelFW_t;
	typedef Likelihood::StateType::Vector::EigenStateOperations operations_t;

	typedef boost::numeric::odeint::euler< stateType_t ,
								   double ,
								   stateType_t,
								   double ,
								   boost::numeric::odeint::vector_space_algebra,
								   operations_t,
								   boost::numeric::odeint::always_resizer > euler_stepper_t;

	static const double N_STEPS_PER_TREE, PROB_THRESHOLD, PROB_TOLERANCE;
	static const size_t REJECTION_SAMPLING_TRY_LIMIT;

	const double FIXED_EULER_DT;

	stateType_t probState;

	kernelsBW_t kernelsBW;
	intKernelBW_t intKernelBW;
	intKernelFW_t intKernelFW;

	Likelihood::Integrator::Base<stateType_t, intKernelBW_t, operations_t>* ptrIntegrator;

	stochasticMappingAlgo_t stochasticMappingAlgo;

	typedef size_t mapStatesKey_t;
	typedef std::pair<stateType_t, stateType_t> mapStatesVal_t;
	typedef std::map< mapStatesKey_t, mapStatesVal_t > mapStates_t;

	std::vector< Likelihood::Scheduler::DAG::NodeDAG* > tasksStack;
	mapStates_t checkpoints;

	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::mean > > accAdaptiveDT;

	void doPreProcessingSteps();
	void doIntegrationStep(Likelihood::Scheduler::DAG::NodeDAG* task);
	void doEventStep(Likelihood::Scheduler::DAG::NodeDAG* task) ;
	void doPostProcessingSteps();
	void doReportState(Likelihood::Scheduler::DAG::NodeDAG* task);

	void doSimulateEvent(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);
	void doSimulateAsyncSpeciationEvent(size_t parentState, Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);

	void doSimulateBranch(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);
	bool doDrawAncestralState(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);
	bool doRejectionSamplingAlgo(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);
	bool doDenseEulerAlgo(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);
	bool doDenseDopriAlgo(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories);

	bool tryDrawingASample(double startTimeFW, double endTimeFW, size_t startStateFW, size_t endStateFW, mapHistories_t::iterator &itHistory);

};

} /* namespace StochasticMapping */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_SEQUENTIALSTOCHASTICMAPPINGCPU_H_ */
