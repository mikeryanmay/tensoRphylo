/*
 * EigenKernels.h
 *
 *  Created on: Sep 2, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_BRANCHWISE_EIGENKERNELS_H_
#define LIKELIHOOD_KERNELS_CPU_BRANCHWISE_EIGENKERNELS_H_

#include "Likelihood/Scheduler/DAG/IncFwdDAG.h"
#include "Likelihood/ConditionTypes/ConditionType.h"

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Data/Reader/IncFwdPhyloReader.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/Approximator/IncFwdLikelihoodApproximator.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Branchwise {

class EigenKernels {
public:
	EigenKernels(Conditions::conditionalProbability_t aConditionType,
				 Phylogeny::Data::ContainerSharedPtr aPtrData,
				 Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager,
				 SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
				 Tensor::ContainerSharedPtr aPtrTensorsContainer,
				 bool &aCondCompatibility);
	~EigenKernels();

	void setInitialCondition(Likelihood::Scheduler::DAG::NodeDAG *node);
	void setInitalExtinctNodeCondition(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void computeAsynchSampling(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void computeAsynchSpeciation(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void computeMassSpeciation(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void computeMassExtinction(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void computeMassSamplingEvent(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void computeMassDestrSamplingEvent(double t, Likelihood::Scheduler::DAG::NodeDAG *node);
	void rescaleProbabilities(Likelihood::Scheduler::DAG::NodeDAG *node);
	double computeLogLikelihood(double t, Eigen::VectorXd &prior, Likelihood::Scheduler::DAG::NodeDAG *node);

private:

	bool &condCompatibility;
	Likelihood::Conditions::conditionalProbability_t conditionType;
	Phylogeny::Data::ContainerSharedPtr ptrData;
	Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr ptrSApproxManager;
	SynchronousEvents::ContainerSharedPtr ptrSynchEventContainer;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

};

} /* namespace Branchwise */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_BRANCHWISE_EIGENKERNELS_H_ */
