/*
 * EigenIntegrationKernelForward.cpp
 *
 *  Created on: May 10, 2020
 *      Author: xaviermeyer
 */

#include "../EigenUtils.h"
#include "../EigenUtilsMatrix.h"
#include "../Misc/StiffnessRatio.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/Approximator/IncLikelihoodApproximator.h"

#include <Eigen/Core>
#include "EigenIntegrationKernelForward.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Specialized {

EigenIntegrationKernelForward::EigenIntegrationKernelForward(
			Likelihood::Approximator::Specialized::SpecializedApproxManagerSharedPtr aPtrSApproxManager,
			Tensor::ContainerSharedPtr aPtrTensorsContainer) :
					ptrSApproxManager(aPtrSApproxManager),
					ptrTensorsContainer(aPtrTensorsContainer) {


}

EigenIntegrationKernelForward::~EigenIntegrationKernelForward() {
}

void EigenIntegrationKernelForward::operator() ( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t ) {
	doIntegrationStep(x, dxdt, t);
}

void EigenIntegrationKernelForward::doIntegrationStep( const Likelihood::StateType::Vector::EigenState &x , Likelihood::StateType::Vector::EigenState &dxdt , double t ) {

	// Recover initial vectors
	Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);
	Eigen::VectorXd vecSum = lambda.col(0) + mu.col(0) + phi.col(0) + delta.col(0);

	// Get U or u hat
	Eigen::VectorXd u;
	ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(t, u);

	// Unobserved: Could be included in main openmp loop with observed x's
	Eigen::VectorXd &dpdt = dxdt.getStateProb();
	const Eigen::VectorXd &p = x.getStateProb();

	// dudt: step 1
	dpdt = -vecSum.cwiseProduct(p);

	// dudt: step 2
	if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_DENSE) {
		Tensor::tensor_t &eta = ptrTensorsContainer->getEta()->getTensor(t);
		dpdt += p.transpose()*eta[0];
	} else if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_SPARSE) {
		Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
		dpdt += p.transpose()*eta[0];
	} else if(ptrTensorsContainer->getEtaStructureType() == Tensor::ETA_QUASSE) {
		Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
		dpdt += defaultComputeQuasseEtaVector(eta(0,1), p);
	}

	// dudt: step 3
	if(ptrTensorsContainer->getOmega()->getDimensions()[0] == 0) { // no tensor
		Eigen::VectorXd resContraction = u.cwiseProduct(p);
		dpdt += 2.*resContraction.cwiseProduct(lambda);
	} else if(ptrTensorsContainer->getOmega()->isSparse()) {
		Eigen::VectorXd resContraction = computeTensorContractionTransposeVectorOpti<true>(ptrTensorsContainer->getOmega(), p, u, t);
		dpdt += resContraction.cwiseProduct(2.*lambda);
	} else {
		assert(false && "Dense omega not allowed.");
	}

	/*
	 * TODO Document that clearly
	 * WARNING This kernel is tailored to work with boost::odeint and decreasing time
	 * We therefore need to invert the sign of the derivative to accomodate for having negative dt
	 *
	 */
	dpdt *= -1.0;

	//std::cout << "dpdt = " << dpdt.transpose() << std::endl;
}

} /* namespace Specialized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
