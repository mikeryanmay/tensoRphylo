/*
 * EigenKernels.h
 *
 *  Created on: Nov 19, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENKERNELS_H_
#define LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENKERNELS_H_

#include "Utils/Parallel/Manager.h"

#if defined(_OPENMP)

#include "Likelihood/StateTypes/OpenMP/EigenState.h"
#include "Likelihood/ConditionTypes/ConditionType.h"

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Data/Reader/IncFwdPhyloReader.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace OpenMP {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, bool withCladoEvents >
class EigenKernels {
public:
	EigenKernels(Phylogeny::Data::ContainerSharedPtr aPtrData,
				 SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
				 Tensor::ContainerSharedPtr aPtrTensorsContainer,
				 bool &aCondCompatibility);
	~EigenKernels();

	void setInitialCondition(const std::vector<PS::Node*>& extantNodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void setInitalExtinctNodeCondition(double t, PS::Node* extinctNode, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void computeAsynchSampling(double t, PS::Node* samplingNode, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void computeAsynchSpeciation(double t, PS::Node* parentNode, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void computeMassSpeciation(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void computeMassExtinction(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void computeMassSamplingEvent(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void computeMassDestrSamplingEvent(double t, const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	void rescaleRequestedBranchesProbabilities(const std::vector<PS::Node*>& nodes, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);
	double computeLogLikelihood(double t, PS::Node* originNode, Eigen::VectorXd &prior, Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x);

private:

	bool &condCompatibility;
	Phylogeny::Data::ContainerSharedPtr ptrData;
	SynchronousEvents::ContainerSharedPtr ptrSynchEventContainer;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

};

} /* namespace OpenMP */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif //defined(_OPENMP)
#endif /* LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENKERNELS_H_ */
