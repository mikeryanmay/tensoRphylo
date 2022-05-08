/*
 * EigenUtils.h
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_EIGENUTILS_H_
#define LIKELIHOOD_KERNELS_CPU_EIGENUTILS_H_

#include <Eigen/Core>

#include "Tensor/IncFwdTensor.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {

template <class T1, class T2>
static inline void computeFirstStepTensorContraction(const T1 &tensor, const T2 &u, Eigen::MatrixXd &resFirstContractionU) {

	size_t dim1 = tensor.size();
	if(dim1 == 0) { // We have nothing: there is no state-dependence and no speciation state change
		resFirstContractionU = u;
	} else { // We have a tensor: there is state-dependence
		size_t dim2 = tensor.front().rows();
		size_t dim3 = tensor.front().cols();
		assert(dim1 == dim2 && dim2 == dim3 && "Tensor must have equal size along each dimensions.");
		for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
			resFirstContractionU.row(iS) = tensor[iS] * u;
		}
	}
}

template <class T1>
static inline void computeSecondStepTensorContractionVector(const Eigen::MatrixXd &omegaTrasposedTimesU, const T1 &vectorP, Eigen::MatrixXd &resContraction) {

	if(omegaTrasposedTimesU.cols() == 1) { // We have nothing: there is no state-dependence and no speciation state change
		resContraction = vectorP.cwiseProduct(omegaTrasposedTimesU.col(0));
	} else { // We have a tensor: there is state-dependence
		resContraction = omegaTrasposedTimesU * vectorP;
	}
}

template <class T1>
static inline void computeSecondStepTensorContractionMatrix(const Eigen::MatrixXd &omegaTrasposedTimesU, const T1 &matrixP, Eigen::MatrixXd &resContraction) {

	if(omegaTrasposedTimesU.cols() == 1) { // We have nothing: there is no state-dependence and no speciation state change
		resContraction = (matrixP.array().colwise() * omegaTrasposedTimesU.col(0).array()).matrix();
	} else { // We have a tensor: there is state-dependence
		resContraction = omegaTrasposedTimesU * matrixP;
	}
}

template <class T1, class T2, class T3>
static inline Eigen::VectorXd computeTensorContractionVector(const T1 &tensor, const T2 &v1, const T3 &v2) {

	size_t dim1 = tensor.size();
	if(dim1 == 0) { // We have nothing: there is no state-dependence and no speciation state change
		Eigen::VectorXd resContraction = v1.cwiseProduct(v2);
		return resContraction;
	} else { // We have a tensor: there is state-dependence
		size_t dim2 = tensor.front().rows();
		size_t dim3 = tensor.front().cols();
		assert(dim1 == dim2 && dim2 == dim3 && "Tensor must have equal size along each dimensions.");
		Eigen::MatrixXd resContraction1(dim2, dim3);
		for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
			//resContraction1.row(iS) = v1.transpose()*tensor[iS];
			resContraction1.row(iS) = tensor[iS]*v2;
		}
		Eigen::VectorXd resContraction2(dim3);
		resContraction2 = resContraction1 * v1;
		return resContraction2;
	}
}

template <bool withCladoEvents, class T1, class T2, class T3>
static inline Eigen::VectorXd computeTensorContractionVectorOpti(const T1 &tensor, const T2 &v1, const T3 &v2, double t) {

	if(!withCladoEvents) { // We have nothing: there is no state-dependence and no speciation state change
		Eigen::VectorXd resContraction = v1.cwiseProduct(v2);
		return resContraction;
	} else { // We have a tensor: there is state-dependence
		const Tensor::sparseTensor_t &omega = tensor->getSparseTensor(t);
		size_t dim1 = omega.size();
		size_t dim2 = omega.front().rows();
		size_t dim3 = omega.front().cols();
		assert(dim1 == dim2 && dim2 == dim3 && "Tensor must have equal size along each dimensions.");
		Eigen::MatrixXd resContraction1(dim2, dim3);
		for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
			//resContraction1.row(iS) = v1.transpose()*omega[iS];
			resContraction1.row(iS) = omega[iS]*v2;
		}
		Eigen::VectorXd resContraction2(dim3);
		resContraction2 = resContraction1 * v1;
		return resContraction2;
	}
}

template <bool withCladoEvents, class T1, class T2, class T3>
static inline Eigen::VectorXd computeTensorContractionTransposeVectorOpti(const T1 &tensor, const T2 &v1, const T3 &v2, double t) {

	if(!withCladoEvents) { // We have nothing: there is no state-dependence and no speciation state change
		Eigen::VectorXd resContraction = v1.cwiseProduct(v2);
		return resContraction;
	} else { // We have a tensor: there is state-dependence
		const Tensor::sparseTensor_t &omega = tensor->getSparseTensor(t);
		size_t dim1 = omega.size();
		size_t dim2 = omega.front().rows();
		size_t dim3 = omega.front().cols();
		assert(dim1 == dim2 && dim2 == dim3 && "Tensor must have equal size along each dimensions.");
		Eigen::MatrixXd resContraction1(dim2, dim3);
		for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
			//resContraction1.row(iS) = v1.transpose()*omega[iS];
			resContraction1.col(iS) = omega[iS]*v2;
		}
		Eigen::VectorXd resContraction2(dim3);
		resContraction2 = resContraction1*v1;
		return resContraction2;
	}
}


} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_EIGENUTILS_H_ */
