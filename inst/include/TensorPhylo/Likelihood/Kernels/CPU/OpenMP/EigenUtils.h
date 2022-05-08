/*
 * EigenUtils.h
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_OPENMP_CPU_EIGENUTILS_H_
#define LIKELIHOOD_KERNELS_OPENMP_CPU_EIGENUTILS_H_

#if defined(_OPENMP)

#include <Eigen/Core>

#include "Tensor/IncFwdTensor.h"

#include "Utils/Parallel/Manager.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace OpenMP {


template <class T1, class T2>
static inline void computeFirstStepTensorContraction(const T1 &tensor, const T2 &u, Eigen::MatrixXd &resFirstContractionU) {

	assert(tensor.size()> 0);

	size_t dim1 = tensor.size();

	if(Utils::Parallel::Manager::getInstance()->useOpenMP()) {
		size_t dim2 = tensor.front().rows();
		size_t dim3 = tensor.front().cols();
		assert(dim1 == dim2 && dim2 == dim3 && "Tensor must have equal size along each dimensions.");

		//omp_set_num_threads(Utils::Parallel::Manager::getInstance()->getMaxNThread());
		#pragma omp parallel
		{
			#pragma omp for schedule(auto)
			for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
				resFirstContractionU.row(iS) = tensor[iS] * u;
			}
		}
	} else {
		size_t dim2 = tensor.front().rows();
		size_t dim3 = tensor.front().cols();
		assert(dim1 == dim2 && dim2 == dim3 && "Tensor must have equal size along each dimensions.");
		for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
			resFirstContractionU.row(iS) = tensor[iS] * u;
		}
	}
}

template <bool withCladoEvents, class T1>
static inline void computeSecondStepTensorContractionVector(const Eigen::MatrixXd &omegaTrasposedTimesU, const T1 &vectorP, Eigen::MatrixXd &resContraction) {

	if(!withCladoEvents) { // We have nothing: there is no state-dependence and no speciation state change
		resContraction = vectorP.cwiseProduct(omegaTrasposedTimesU.col(0));
	} else { // We have a tensor: there is state-dependence
		resContraction = omegaTrasposedTimesU * vectorP;
	}
}

template <class T1>
static inline void computeSecondStepTensorContractionMatrix(const Eigen::MatrixXd &omegaTrasposedTimesU, const T1 &matrixP, Eigen::MatrixXd &resContraction) {

	if(omegaTrasposedTimesU.cols() == 1) { // We have nothing: there is no state-dependence and no speciation state change

		if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(matrixP.cols())) {
			resContraction.resize(matrixP.rows(), matrixP.cols());

			size_t nThreads = Utils::Parallel::Manager::getInstance()->getNThread();
			std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineChunks(nThreads, matrixP.cols()));

			//omp_set_num_threads(Utils::Parallel::Manager::getInstance()->getMaxNThread());
			#pragma omp parallel
			{
				//omp_set_num_threads(nThreads);

				#pragma omp for schedule(static)
				for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
					resContraction.block(0, chunks[iT], matrixP.rows(), chunks[iT+1]-chunks[iT]) = (matrixP.block(0, chunks[iT], matrixP.rows(), chunks[iT+1]-chunks[iT]).array().colwise() * omegaTrasposedTimesU.col(0).array()).matrix();
				}
			}
		} else {
			resContraction = (matrixP.array().colwise() * omegaTrasposedTimesU.col(0).array()).matrix();
		}
	} else { // We have a tensor: there is state-dependence
		if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(matrixP.cols())) {
			resContraction.resize(matrixP.rows(), matrixP.cols());

			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getNThread();
			std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, matrixP.cols()));
			size_t nThreads = chunks.size()-1;

			//omp_set_num_threads(nThreads);
			#pragma omp parallel
			{
				#pragma omp for schedule(static)
				for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
					resContraction.block(0, chunks[iT], matrixP.rows(), chunks[iT+1]-chunks[iT]) = omegaTrasposedTimesU * matrixP.block(0, chunks[iT], matrixP.rows(), chunks[iT+1]-chunks[iT]);
				}
			}
		} else {
			resContraction = omegaTrasposedTimesU * matrixP;
		}
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

		if(Utils::Parallel::Manager::getInstance()->useOpenMP()) {
			//omp_set_num_threads(Utils::Parallel::Manager::getInstance()->getMaxNThread());
			#pragma omp parallel
			{
				#pragma omp for schedule(auto)
				for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
					resContraction1.row(iS) = v1.transpose() * tensor[iS];
				}

			}
		} else {
			for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
				resContraction1.row(iS) = v1.transpose() * tensor[iS];
			}
		}

		Eigen::VectorXd resContraction2(dim3);
		resContraction2 = resContraction1 * v2;
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
		if(Utils::Parallel::Manager::getInstance()->useOpenMP()) {
			#pragma omp parallel
			{
				#pragma omp for schedule(auto)
				for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
					resContraction1.row(iS) = v1.transpose() *  omega[iS];
				}

			}
		} else {
			for(size_t iS=0; iS<dim1; ++iS) { // Potential for pragma openmp
				resContraction1.row(iS) = v1.transpose() * omega[iS];
			}
		}
		Eigen::VectorXd resContraction2(dim3);
		resContraction2 = resContraction1 * v2;
		return resContraction2;
	}
}

} /* namespace OpenMP */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif

#endif /* LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENUTILS_H_ */
