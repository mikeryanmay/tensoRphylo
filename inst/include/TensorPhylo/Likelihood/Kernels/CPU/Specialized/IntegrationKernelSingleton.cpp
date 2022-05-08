/*
 * IntegrationKernelSingleton.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "IntegrationKernelSingleton.h"

#include "../EigenUtils.h"
#include "../EigenUtilsMatrix.h"
#include "Tensor/IncTensor.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/Approximator/Specialized/DenseUnobserved.h"

#include <Eigen/Core>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

IntegrationKernelSingleton::IntegrationKernelSingleton(
			Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
			Likelihood::Approximator::Specialized::DenseUnobservedSharedPtr  aPtrUnobsProb,
			Tensor::ContainerSharedPtr aPtrTensorsContainer) :
					PROCESS_TYPE(aExtantProcessType),
					ptrDenseUnobs(aPtrUnobsProb),
					ptrTensorsContainer(aPtrTensorsContainer) {

	assert(PROCESS_TYPE == ptrDenseUnobs->getProcessType());
	resFirstContractionU.resize(ptrTensorsContainer->getNumberOfState(), ptrTensorsContainer->getNumberOfState());

}

IntegrationKernelSingleton::~IntegrationKernelSingleton() {
}

void IntegrationKernelSingleton::operator() ( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t ) {
	doIntegrationStep(x, dxdt, t);
}

void IntegrationKernelSingleton::doIntegrationStep( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t ) {

	// Recover initial vectors
	Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);
	Eigen::VectorXd vecSum;

	// Get U or u hat
	Eigen::VectorXd u;
	ptrDenseUnobs->defineProbabilities(t, u);

	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		vecSum = lambda.col(0) + mu.col(0) + phi.col(0) + delta.col(0);
	} else {
		vecSum = lambda.col(0) + mu.col(0);
	}

	// Unobserved: Could be included in main openmp loop with observed x's
	Eigen::VectorXd &dsdt = dxdt.getStateProb();
	const Eigen::VectorXd &s = x.getStateProb();

	// dsdt: step 1
	dsdt = - vecSum.cwiseProduct(s);

	// dsdt: step 2
	if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_DENSE) {
		Tensor::tensor_t &eta = ptrTensorsContainer->getEta()->getTensor(t);
		dsdt += eta[0]*s;
	} else if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_SPARSE) {
		Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
		dsdt += eta[0]*s;
	} else if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_QUASSE) {
		Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
		dsdt += defaultComputeQuasseEtaVector(eta(0,1), s);
	}

	// dsdt: step 3
	Eigen::MatrixXd resContraction;
	if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
		resContraction = s.cwiseProduct(u);
		dsdt += 2.0 * resContraction.cwiseProduct(lambda);
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		//Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		//resContraction = computeTensorContractionVector(omega, u, s);
		Eigen::VectorXd resContraction = computeTensorContractionVectorOpti<true>(ptrTensorsContainer->getOmega(), s, u, t);
		dsdt += resContraction.cwiseProduct(2.*lambda);
	} else {
		assert(false && "Dense omega not allowed.");
	}


	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		dsdt += phi.cwiseProduct(u) + delta;
	}


}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
