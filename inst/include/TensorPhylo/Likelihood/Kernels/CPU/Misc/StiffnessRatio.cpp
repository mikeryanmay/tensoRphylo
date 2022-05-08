/*
 * StiffnessRatio.cpp
 *
 *  Created on: Apr 9, 2020
 *      Author: meyerx
 */

#include "StiffnessRatio.h"

#include "../EigenUtils.h"
#include "Tensor/IncTensor.h"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace Likelihood {
namespace Kernels {
namespace CPU {

StiffnessRatio::StiffnessRatio(Tensor::ContainerSharedPtr aPtrTensorsContainer) :
		ptrTensors(aPtrTensorsContainer) {
}

StiffnessRatio::~StiffnessRatio() {
}

std::pair<double, bool> StiffnessRatio::estimateStiffnessRatio(double t, const Eigen::VectorXd &u) {

	// Matrix is eta
	Eigen::MatrixXd A = ptrTensors->getEigenMatrixEta(t);

	// Diagonal elements
	Eigen::VectorXd diag = -(ptrTensors->getEigenVecLambda(t).col(0) + ptrTensors->getEigenVecMu(t).col(0) + ptrTensors->getEigenVecPhi(t).col(0) + ptrTensors->getEigenVecDelta(t).col(0));

	A.diagonal() += diag;

	if(ptrTensors->getOmega()->getDimensions()[0] == 0) {
		A.diagonal() += 2*ptrTensors->getEigenVecLambda(t).cwiseProduct(u);
	} else {

		Eigen::VectorXd uTimesLambda = 2.0*ptrTensors->getEigenVecLambda(t).col(0).cwiseProduct(u);
		Eigen::MatrixXd stepContractionU(u.size(), u.size());
		computeFirstStepTensorContraction(ptrTensors->getOmega()->getSparseTensor(t), uTimesLambda, stepContractionU);
		A += stepContractionU;
	}

	Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
	Eigen::VectorXd reLambda = solver.eigenvalues().real();
	double maxReLambda = reLambda.cwiseAbs().maxCoeff();
	double minReLambda = reLambda.cwiseAbs().minCoeff();

	/*std::cout << "t = " << t << std::endl;
	std::cout << "u : " << std::endl;
	std::cout << u.transpose() << std::endl;
	std::cout << "A:" << std::endl;
	std::cout << A << std::endl;
	std::cout << "eig(A) - real :" << std::endl;
	std::cout << solver.eigenvalues().real().transpose() << std::endl;*/
	//std::cout << "eig(A) - im :" << std::endl;
	//std::cout << solver.eigenvalues().imag().transpose() << std::endl;

	//std::cout << solver.eigenvalues().imag().cwiseAbs().maxCoeff() << std::endl;
	return std::make_pair(maxReLambda/minReLambda, solver.eigenvalues().imag().cwiseAbs().maxCoeff() > 1.e-7);

}

} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
