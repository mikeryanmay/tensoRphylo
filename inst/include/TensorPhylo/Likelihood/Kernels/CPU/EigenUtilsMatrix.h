/*
 * EigenUtilsMatrix.h
 *
 *  Created on: April 1, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_KERNELS_CPU_EIGENUTILSMATRIX_H_
#define LIKELIHOOD_KERNELS_CPU_EIGENUTILSMATRIX_H_

#include <iostream>
#include <Eigen/Core>

#include "Tensor/IncFwdTensor.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {


// Tentative implementation of a fast "convolution 1D" for QUASSE type matrix
// These matrices are equivalent to filters of form [alpha, -2*alpha, alpha]
template <int K, class T1>
static Eigen::VectorXd computeQuasseEtaVector(double alpha, const T1 &v1) {

	const size_t N = v1.size();
	size_t M = std::floor((float)N/(float)K);
	M = (M*K == N) ? M-1 : M;
	Eigen::VectorXd dataOut(N);
	Eigen::Matrix<double, K+2, 1> buffer;

	assert(N>= 3*1 && "Use a lower stride (K) or a sparse/dense eigen::matrix (small state-space problem).");

	{ 	// Border
		// 0 to K
		//Eigen::Matrix<double, K+1, 1> buffer;
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
		for(size_t j=0; j<K+1; j++) {
			buffer(j) = alpha*v1(j);
		}

		// compute y = a*x(-1) -2*a*x() + a*x(+1)
		dataOut(0) = -buffer[0]+buffer[1];
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
		for(size_t j=1; j<K; j++) {
			dataOut(j) = buffer(j-1);
			dataOut(j) -= 2.*buffer(j);
			dataOut(j) += buffer(j+1);
		}
	}

	{	// Inner computation
		// K to M*K
		for(size_t i=K; i<M*K; i+=K) {
			 // preload alpha*x
			//Eigen::Matrix<double, K+2, 1> buffer;
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
			for(size_t j=0; j<K+2; j++) {
				buffer(j) = alpha*v1(i+j-1);
			}

			// compute y = a*x(-1) -2*a*x() + a*x(+1)
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
			for(size_t j=0; j<K; j++) {
				dataOut(i+j) = buffer(j);
				dataOut(i+j) -= 2.*buffer(j+1);
				dataOut(i+j) += buffer(j+2);
			}
		}
	}

	{	 // Border
		// M*K to N
		const size_t L=N-M*K;
		//Eigen::MatrixXd buffer(L+1, 1);
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
		for(size_t j=0; j<L+1; j++) {
			buffer(j) = alpha*v1(M*K+j-1);
		}

		// compute y = a*x(-1) -2*a*x() + a*x(+1)
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
		for(size_t j=0; j<L-1; j++) {
			dataOut(M*K+j) = buffer(j);
			dataOut(M*K+j) -= 2.*buffer(j+1);
			dataOut(M*K+j) += buffer(j+2);
		}
		dataOut(N-1) = buffer(L-1)-buffer(L);
	}

	return dataOut;
}

template <class T1>
static inline Eigen::VectorXd defaultComputeQuasseEtaVector(double alpha, const T1 &v1) {
	return computeQuasseEtaVector< 4 >(alpha, v1);
}


// Tentative implementation of a fast "convolution 1D" for QUASSE type matrix
// These matrices are equivalent to filters of form [alpha, -2*alpha, alpha]
// Here we extend that to simultaneous application to multiple vectors
// TODO The outer loop would be a perfect candidate for openMP
template <int K1, int K2, class T1>
static Eigen::MatrixXd computeQuasseEtaMatrix(double alpha, const T1 &m1) {

	const size_t N_R = m1.rows();
	size_t M_R = std::floor((float)N_R/(float)K1);
	M_R = M_R*K1 == N_R ? M_R-1 : M_R;

	assert(N_R>= 3*K1 && "Use a lower stride (K1) or a sparse/dense eigen::matrix (small state-space problem).");

	const size_t N_C = m1.cols();
	const size_t M_C = std::floor((float)N_C/(float)K2);

	Eigen::MatrixXd dataOut(N_R, N_C);
	Eigen::Matrix<double, K1+2, K2> buffer;

	for(size_t iC=0; iC<M_C*K2; iC+=K2) {
		{	// Border
			// 0 to K1
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
			for(size_t jC=0; jC<K2; ++jC) {
				//#pragma GCC ivdep
				for(size_t jR=0; jR<K1+1; jR++) {
					buffer(jR, jC) = alpha*m1(jR, iC+jC);
				}

				// compute y = a*x(-1) -2*a*x() + a*x(+1)
				dataOut(0, iC+jC) = -buffer(0, jC)+buffer(1, jC);
				//#pragma GCC ivdep
				for(size_t jR=1; jR<K1; ++jR) {
					dataOut(jR, iC+jC) = buffer(jR-1, jC);
					dataOut(jR, iC+jC) -= 2.*buffer(jR, jC);
					dataOut(jR, iC+jC) += buffer(jR+1, jC);
				}
			}
		}


		{	// Inner computation
			// K1 to M*K1
			for(size_t iR=K1; iR<M_R*K1; iR+=K1) {
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
				for(size_t jC=0; jC<K2; ++jC) {
					//#pragma GCC ivdep
					for(size_t jR=0; jR<K1+2; jR++) {
						buffer(jR, jC) = alpha*m1(iR+jR-1, iC+jC);
					}

					// compute y = a*x(-1) -2*a*x() + a*x(+1)
					//#pragma GCC ivdep
					for(size_t jR=0; jR<K1; ++jR) {
						dataOut(iR+jR, iC+jC) = buffer(jR, jC);
						dataOut(iR+jR, iC+jC) -= 2.*buffer(jR+1, jC);
						dataOut(iR+jR, iC+jC) += buffer(jR+2, jC);
					}
				}
			}
		}

		{	// Border
			// M_R*K1 to N_R
			const size_t L_R = N_R-M_R*K1;
#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
			for(size_t jC=0; jC<K2; ++jC) {
				//#pragma GCC ivdep
				for(size_t jR=0; jR<L_R+1; jR++) {
					buffer(jR, jC) = alpha*m1(M_R*K1+jR-1, iC+jC);
				}

				// compute y = a*x(-1) -2*a*x() + a*x(+1)
				//#pragma GCC ivdep
				for(size_t jR=0; jR<L_R-1; ++jR) {
					dataOut(M_R*K1+jR, iC+jC) = buffer(jR, jC);
					dataOut(M_R*K1+jR, iC+jC) -= 2.*buffer(jR+1, jC);
					dataOut(M_R*K1+jR, iC+jC) += buffer(jR+2, jC);
				}
				dataOut(N_R-1, iC+jC) = buffer(L_R-1, jC)-buffer(L_R, jC);
			}
		}
	}

	const size_t L_C = N_C-M_C*K2;
	if(L_C > 0) {		// Last columns
		{	// Border
			// 0 to K1
#if (defined(__GNUC__) && !defined(__clang__))
			#pragma GCC ivdep
#endif
			for(size_t jC=0; jC<L_C; ++jC) {
				//#pragma GCC ivdep
				for(size_t jR=0; jR<K1+1; jR++) {
					buffer(jR, jC) = alpha*m1(jR, M_C*K2+jC);
				}

				// compute y = a*x(-1) -2*a*x() + a*x(+1)
				dataOut(0, M_C*K2+jC) = -buffer(0, jC)+buffer(1, jC);
				//#pragma GCC ivdep
				for(size_t jR=1; jR<K1; ++jR) {
					dataOut(jR, M_C*K2+jC) = buffer(jR-1, jC);
					dataOut(jR, M_C*K2+jC) -= 2.*buffer(jR, jC);
					dataOut(jR, M_C*K2+jC) += buffer(jR+1, jC);
				}
			}
		}


		{	// Inner computation
			// K1 to M*K1
			for(size_t iR=K1; iR<M_R*K1; iR+=K1) {

#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
				for(size_t jC=0; jC<L_C; ++jC) {
					//#pragma GCC ivdep
					for(size_t jR=0; jR<K1+2; jR++) {
						buffer(jR, jC) = alpha*m1(iR+jR-1, M_C*K2+jC);
					}

					// compute y = a*x(-1) -2*a*x() + a*x(+1)
					//#pragma GCC ivdep
					for(size_t jR=0; jR<K1; ++jR) {
						dataOut(iR+jR, M_C*K2+jC) = buffer(jR, jC);
						dataOut(iR+jR, M_C*K2+jC) -= 2.*buffer(jR+1, jC);
						dataOut(iR+jR, M_C*K2+jC) += buffer(jR+2, jC);
					}
				}
			}
		}

		{	// Border
			// M_R*K1 to N_R
			const size_t L_R = N_R-M_R*K1;

#if (defined(__GNUC__) && !defined(__clang__))
		#pragma GCC ivdep
#endif
			for(size_t jC=0; jC<L_C; ++jC) {
				//#pragma GCC ivdep
				for(size_t jR=0; jR<L_R+1; jR++) {
					buffer(jR, jC) = alpha*m1(M_R*K1+jR-1, M_C*K2+jC);
				}

				// compute y = a*x(-1) -2*a*x() + a*x(+1)
				//#pragma GCC ivdep
				for(size_t jR=0; jR<L_R-1; ++jR) {
					dataOut(M_R*K1+jR, M_C*K2+jC) = buffer(jR, jC);
					dataOut(M_R*K1+jR, M_C*K2+jC) -= 2.*buffer(jR+1, jC);
					dataOut(M_R*K1+jR, M_C*K2+jC) += buffer(jR+2, jC);
				}
				dataOut(N_R-1, M_C*K2+jC) = buffer(L_R-1, jC)-buffer(L_R, jC);
			}
		}
	}

	return dataOut;
}


template <class T1>
static inline Eigen::MatrixXd defaultComputeQuasseEtaMatrix(double alpha, const T1 &m1) {
	return computeQuasseEtaMatrix< 8, 8 >(alpha, m1);
}

} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_EIGENUTILSMATRIX_H_ */
