/*
 * IntegrationKernelPairedUS.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "IntegrationKernelPairedUS.h"

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

IntegrationKernelPairedUS::IntegrationKernelPairedUS(
			Likelihood::Approximator::Specialized::processType_t aExtantProcessType,
			Tensor::ContainerSharedPtr aPtrTensorsContainer) :
					PROCESS_TYPE(aExtantProcessType),
					ptrTensorsContainer(aPtrTensorsContainer) {

	resFirstContractionU.resize(ptrTensorsContainer->getNumberOfState(), ptrTensorsContainer->getNumberOfState());

}

IntegrationKernelPairedUS::~IntegrationKernelPairedUS() {
}

void IntegrationKernelPairedUS::operator() ( const Likelihood::StateType::Matrix::EigenState &x , Likelihood::StateType::Matrix::EigenState &dxdt , double t ) {
	doIntegrationStep(x, dxdt, t);
}

void IntegrationKernelPairedUS::doIntegrationStep( const Likelihood::StateType::Matrix::EigenState &x , Likelihood::StateType::Matrix::EigenState &dxdt , double t ) {

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

	// Getting prob vectors
	Eigen::VectorXd &dudt = dxdt.getStateProb(0);
	const Eigen::VectorXd &u = x.getStateProb(0);
	Eigen::VectorXd &dsdt = dxdt.getStateProb(1);
	const Eigen::VectorXd &s = x.getStateProb(1);

	// singleton
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
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		resContraction = computeTensorContractionVector(omega, u, s);
	} else {
		assert(false && "Dense omega not allowed.");
	}
	dsdt += 2.0 * resContraction.cwiseProduct(lambda);


	if(PROCESS_TYPE == Likelihood::Approximator::Specialized::FULL_PROCESS) {
		dsdt += phi.cwiseProduct(u) + delta;
	}


	// Unobserved

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
		Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		resContraction = computeTensorContractionVector(omega, u, u);
		dudt += resContraction.cwiseProduct(lambda);
	} else {
		assert(false && "Dense omega not allowed.");
	}



}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
