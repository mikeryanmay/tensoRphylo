/*
 * QuasistationaryFrequency.cpp
 *
 *  Created on: May 6, 2022
 *      Author: mrmay
 */

#include <Eigen/Dense>
#include <limits>
#include "Likelihood/Kernels/CPU/Misc/QuasistationaryFrequency.h"
#include "Tensor/IncTensor.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {

Eigen::VectorXd computeQuasiStationaryFrequency(Tensor::ContainerSharedPtr ptrTensors, double t) {

	// get diversification parameters
	Eigen::MatrixXd &mu     = ptrTensors->getEigenVecMu(t);
	Eigen::MatrixXd &lambda = ptrTensors->getEigenVecLambda(t);

	// create the amplification matrix
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(mu.size(), mu.size());

	// subtract extinction rate
	A.diagonal() = -1.0 * mu;

	// apply anagenetic events
	if (ptrTensors->getEtaStructureType() == Tensor::ETA_DENSE) {
		Tensor::tensor_t &eta = ptrTensors->getEta()->getTensor(t);
	  A += eta[0];
	} else if (ptrTensors->getEtaStructureType() == Tensor::ETA_SPARSE) {
		Tensor::sparseTensor_t &eta = ptrTensors->getEta()->getSparseTensor(t);
	  A += eta[0];
	} else if (ptrTensors->getEtaStructureType() == Tensor::ETA_QUASSE) {
		Tensor::tensor_t &eta = ptrTensors->getEta()->getTensor(t);
	  A += eta[0];
	}

	// apply cladogenetic events
	if (ptrTensors->getOmega()->getDimensions()[0] == 0) { // no tensor
		// just add lambda to diagonal elements
		A.diagonal() += lambda;
	} else if (ptrTensors->getOmega()->isSparse()) {

		// create the lambda component
		Eigen::MatrixXd L = Eigen::MatrixXd::Zero( lambda.size(), lambda.size() );

		// get omega for this time point
		Tensor::sparseTensor_t &omega = ptrTensors->getOmega()->getSparseTensor(t);

		// loop over sparse matrices in omega
		for(size_t i = 0; i < omega.size(); ++i) {

			// get this lambda value
			double this_lambda = lambda(i,0);

			// get the sparse matrix
			Tensor::eigenSparseMatrix_t &slice = omega.at(i);

			// loop over rows
			for (size_t k = 0; k < slice.outerSize(); ++k) {

				// loop over elements of the row (columns)
				for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(slice, k); it; ++it) {

					// increment appropriate Lambda values
					double val = it.value();
					size_t row = it.row();
					size_t col = it.col();
					if ( row == col && row == i ) {

						// birth into its own state
						L(row, i) += this_lambda * val;

						// birth events not incrementing this state
						L(row, i) -= this_lambda * (1.0 - val);

					} else {

						// birth from one state to another
						L(row, i) += this_lambda * val;
						L(col, i) += this_lambda * val;

					}

			  } // end loop over elements of this row (left daughters)
			} // end loop over rows (right daughters)
		} // end loop over ancestor states

		// add Lambda to the amplification matrix
		A += L;

	} else {
		assert(false && "Dense omega not allowed.");
	}

	// do eigen decomposition
	Eigen::EigenSolver<Eigen::MatrixXd> eigen;
	eigen.compute(A);
	Eigen::VectorXd eigen_values  = eigen.eigenvalues().real();
	Eigen::MatrixXd eigen_vectors = eigen.eigenvectors().real();

	// find the largest eigen value
	double max = -std::numeric_limits<double>::max();
	size_t largest_index;
	for (size_t i = 0; i < lambda.size(); ++i) {
		if ( eigen_values(i) > max ) {
			max = eigen_values(i);
			largest_index = i;
		}
	}

	// get the corresponding eigenvector and normalize
	Eigen::VectorXd dominant_eigenvector = eigen_vectors.col(largest_index);
	dominant_eigenvector /= dominant_eigenvector.sum();

	return dominant_eigenvector;

}


} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */
