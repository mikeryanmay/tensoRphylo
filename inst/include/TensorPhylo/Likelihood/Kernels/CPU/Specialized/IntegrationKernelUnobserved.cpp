/*
 * IntegrationKernelUnobserved.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "IntegrationKernelUnobserved.h"

#include "../EigenUtils.h"
#include "../EigenUtilsMatrix.h"
#include "../Misc/StiffnessRatio.h"
#include "Tensor/IncTensor.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"

#include <Eigen/Core>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

IntegrationKernelUnobserved::IntegrationKernelUnobserved(
		Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
		Tensor::ContainerSharedPtr aPtrTensorsContainer) :
				PROCESS_TYPE(aExtantProcessType), ptrTensorsContainer(aPtrTensorsContainer) {

	resFirstContractionU.resize(ptrTensorsContainer->getNumberOfState(), ptrTensorsContainer->getNumberOfState());

}

IntegrationKernelUnobserved::~IntegrationKernelUnobserved() {
}

void IntegrationKernelUnobserved::operator() ( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t ) {
	doIntegrationStep(x, dxdt, t);
}

void IntegrationKernelUnobserved::doIntegrationStep( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t ) {

	// Recover initial vectors
	Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);
	Eigen::VectorXd vecSum;

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		vecSum = lambda.col(0) + mu.col(0) + phi.col(0) + delta.col(0);
	} else {
		vecSum = lambda.col(0) + mu.col(0);
	}

	// Unobserved: Could be included in main openmp loop with observed x's
	Eigen::VectorXd &dudt = dxdt.getStateProb();
	const Eigen::VectorXd &u = x.getStateProb();

	// dudt: step 1
	dudt = mu - vecSum.cwiseProduct(u);

	// dudt: step 2
	if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_DENSE) {
		Eigen::MatrixXd eta = ptrTensorsContainer->getEigenMatrixEta(t);
		dudt += eta*u;
	} else if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_SPARSE) {
		Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
		dudt += eta[0]*u;
	} else if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_QUASSE) {
		Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
		dudt += defaultComputeQuasseEtaVector(eta(0,1), u);
	}

	// dudt: step 3
	if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
		Eigen::VectorXd resContraction = u.cwiseProduct(u);
		dudt += resContraction.cwiseProduct(lambda);
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		Eigen::VectorXd resContraction = computeTensorContractionVectorOpti<true>(ptrTensorsContainer->getOmega(), u, u, t);
		dudt += resContraction.cwiseProduct(lambda);
	} else {
		assert(false && "Dense omega not allowed.");
	}

}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
