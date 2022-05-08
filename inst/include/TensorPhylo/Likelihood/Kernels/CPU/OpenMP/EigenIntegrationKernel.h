/*
 * EigenInterationKernel.h
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENINTEGRATIONKERNEL_H_
#define LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENINTEGRATIONKERNEL_H_

#include "Utils/Parallel/Manager.h"
#if defined(_OPENMP)

#include "Utils/Profiling/CustomProfiling.h"
#include "Utils/Profiling/UniqueProfiler.h"
#include "Likelihood/StateTypes/OpenMP/EigenState.hpp"

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Data/Reader/IncFwdPhyloReader.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/ConditionTypes/ConditionType.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace OpenMP {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
class EigenIntegrationKernel {
public:
	EigenIntegrationKernel(const size_t N_MAX_STATE_VECTOR,
							Phylogeny::Data::ContainerSharedPtr aPtrData,
							SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
							Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~EigenIntegrationKernel();

	void operator() ( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t );

	std::pair<double, bool> getStiffnessRatio(double t, const Eigen::VectorXd &u) const;

private:

	bool isPrecomputedEtaAvailable;
	Tensor::eigenSparseMatrix_t precomputedEtaSparse;
	Eigen::MatrixXd resFirstContractionU, precomputedEta;

	Phylogeny::Data::ContainerSharedPtr ptrData;
	SynchronousEvents::ContainerSharedPtr ptrSynchEventContainer;
	Tensor::ContainerSharedPtr ptrTensorsContainer;

	//void doIntegrationStepChunk( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t );
	void doIntegrationStepBlock( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t );
	//void doIntegrationStepChunk2GEMM( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t );
	void doIntegrationStepBlock2GEMM( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t );

	void computeConditioningStateProb( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t );
};

} /* namespace OpenMP */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif //defined(_OPENMP)
#endif /* LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENINTEGRATIONKERNEL_H_ */
