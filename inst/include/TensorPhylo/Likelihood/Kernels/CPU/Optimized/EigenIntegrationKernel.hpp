/*
 * EigenIntegrationKernel.cpp
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#include "EigenIntegrationKernel.h"

#include "../EigenUtils.h"
#include "../EigenUtilsMatrix.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../Misc/StiffnessRatio.h"
#include "Utils/Profiling/UniqueProfiler.h"

#ifndef LIKELIHOOD_KERNELS_CPU_OPTIMIZED_EIGENINTEGRATIONKERNEL_DEF_H_
#define LIKELIHOOD_KERNELS_CPU_OPTIMIZED_EIGENINTEGRATIONKERNEL_DEF_H_

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace Optimized {

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::EigenIntegrationKernel(
			const size_t N_MAX_STATE_VECTOR,
			Phylogeny::Data::ContainerSharedPtr aPtrData,
			SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
			Tensor::ContainerSharedPtr aPtrTensorsContainer) :
					ptrData(aPtrData),
					ptrSynchEventContainer(aPtrSynchEventContainer),
					ptrTensorsContainer(aPtrTensorsContainer) {

	// PROFILING tStep = 0.;
	resFirstContractionU.resize(ptrTensorsContainer->getNumberOfState(), ptrTensorsContainer->getNumberOfState());

	isPrecomputedEtaAvailable = ptrTensorsContainer->getEta()->isContantThroughTime() &&
								ptrTensorsContainer->getMu()->isContantThroughTime() &&
								ptrTensorsContainer->getLambda()->isContantThroughTime() &&
								ptrTensorsContainer->getDelta()->isContantThroughTime() &&
								ptrTensorsContainer->getPhi()->isContantThroughTime() &&
								etaStructure != Tensor::ETA_QUASSE;

	if(isPrecomputedEtaAvailable) {
		//std::cout << "Use precomputed eta" << std::endl;
		double t=0.;
		const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
		const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
		const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
		const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);
		if(etaStructure == Tensor::ETA_SPARSE) {
			precomputedEtaSparse = ptrTensorsContainer->getEta()->getSparseTensor(t)[0];
			precomputedEtaSparse.diagonal() -= mu + lambda + phi + delta;
			precomputedEtaSparse.makeCompressed();
		} else if(etaStructure == Tensor::ETA_DENSE) {
			precomputedEta = ptrTensorsContainer->getEta()->getTensor(t)[0];
			precomputedEta.diagonal() -= mu + lambda + phi + delta;
		}
	}
}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::~EigenIntegrationKernel() {
}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::operator() ( const Likelihood::StateType::Optimized::EigenState<conditionalProbType> &x , Likelihood::StateType::Optimized::EigenState<conditionalProbType> &dxdt , double t ) {
	doIntegrationStep(x, dxdt, t);

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
std::pair<double, bool> EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::getStiffnessRatio(double t, const Eigen::VectorXd &u) const {
	StiffnessRatio sr(ptrTensorsContainer);
	return sr.estimateStiffnessRatio(t, u);
}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::doIntegrationStep( const Likelihood::StateType::Optimized::EigenState<conditionalProbType> &x , Likelihood::StateType::Optimized::EigenState<conditionalProbType> &dxdt , double t ) {

	// Recover initial vectors
	const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	const Eigen::Ref< const Eigen::VectorXd > u = x.getUnobservedStateProb();

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1")
	// For all probability vector: - (mu+delta+mu+lambda) * [u/p/...]+ eta*[u/p/...]
	{
		Eigen::Ref< Eigen::MatrixXd > dAlldt = dxdt.getStateProb();
		const Eigen::Ref< const Eigen::MatrixXd > all = x.getStateProb();

		if(etaStructure == Tensor::ETA_SPARSE) {
			if(isPrecomputedEtaAvailable) {
				dAlldt = precomputedEtaSparse*all;
			} else {
				const Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
				Eigen::MatrixXd vecSum = phi + delta + mu + lambda;
				dAlldt = eta[0]*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
			}
		} else if(etaStructure == Tensor::ETA_DENSE) {
			if(isPrecomputedEtaAvailable) {
				dAlldt = precomputedEta*all;
			} else {
				const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
				Eigen::MatrixXd vecSum = phi + delta + mu + lambda;
				dAlldt = eta*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
			}
		} else if(etaStructure == Tensor::ETA_QUASSE) {
			const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
			double scaledSigma = eta(0,1);
			Eigen::MatrixXd vecSum = phi + delta + mu + lambda;
			dAlldt = defaultComputeQuasseEtaMatrix(scaledSigma, all) - (all.array().colwise()* vecSum.col(0).array()).matrix();
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1")

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")
	// For all [u/p]: d?dt += cst * lambda * omega * [p/u] u
	{
		Eigen::Ref< Eigen::MatrixXd > dupdt = dxdt.getUnobservedAndObservedStateProb();
		const Eigen::Ref< const Eigen::MatrixXd > up = x.getUnobservedAndObservedStateProb();
		Eigen::MatrixXd resContractionMatrix(up.rows(), up.cols());
		Eigen::ArrayXd tmpArray1 = lambda.col(0).array();
		Eigen::ArrayXd tmpArray2 = 2.*Eigen::ArrayXd::Ones(dupdt.cols());
		tmpArray2(0) = 1.;

		if(withCladoEvents) {
			const Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
			/*for(size_t iS=0; iS<omega.size(); ++iS) { // Potential for pragma openmp
				resFirstContractionU.col(iS) = omega[iS] * u;
			}
			dupdt += (((resFirstContractionU * up).array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();*/
			computeFirstStepTensorContraction(omega, u, resFirstContractionU);
			//std::cout << "Matrix resContr at t = " <<  t << " : " << std::endl << resFirstContractionU << std::endl << "--------------------------------------" << std::endl;
			resContractionMatrix = resFirstContractionU * up;
			dupdt += ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
		} else {
			resFirstContractionU = u;
			resContractionMatrix = (up.array().colwise() * resFirstContractionU.col(0).array()).matrix();
			dupdt += ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
			//dupdt += ((up.array().colwise() * (u.array()*tmpArray1)).rowwise() * tmpArray2.transpose()).matrix();
		}
		_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")

	}

	// Unobserved: += mu
	{
		Eigen::Ref< Eigen::VectorXd > dudt = dxdt.getUnobservedStateProb();
		dudt += mu;
	}

	// singleton:
	if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > dsdt = dxdt.getSingletonStateProb();
		const Eigen::Ref< const Eigen::VectorXd > s = x.getSingletonStateProb();
		// dsdt: step 1 - correction for phi since dsdt -= (mu + lambda + delta)*s
		// dsdt: step 1 - but we did dsdt -= (phi + mu + lambda + delta)*s
//		dsdt += phi.cwiseProduct(s);
		// dsdt: step 2
		if(withCladoEvents) {
			dsdt += 2.0 * (resFirstContractionU * s).cwiseProduct(lambda);
		} else {
			dsdt += 2.0 * (u.array() * s.array() * lambda.array()).matrix();
		}

		// dsdt: step 3
		dsdt += (phi).cwiseProduct(u) + delta;
	}

	// unobserved no sampling:
	if(conditionalProbType == Conditions::ROOT_SURVIVAL || conditionalProbType == Conditions::ROOT_MRCA ||
	   conditionalProbType == Conditions::STEM_SURVIVAL) {
		Eigen::Ref< Eigen::VectorXd > dudt_hat = dxdt.getUnobservedNoSamplingStateProb();
		const Eigen::Ref< const Eigen::VectorXd > u_hat = x.getUnobservedNoSamplingStateProb();
		// duhatdt: step 1 - correction for phi since duhatdt -= (mu + lambda)*u_hat
		// duhatdt: step 1 - but we did duhatdt -= (phi + mu + lambda + delta)*u_hat
		// add mu too
		dudt_hat += mu + (phi + delta).cwiseProduct(u_hat);

		// dudt: step 3
		Eigen::VectorXd resContraction = computeTensorContractionVectorOpti<withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< const Eigen::VectorXd >, Eigen::Ref< const Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), u_hat, u_hat, t);
		dudt_hat += resContraction.cwiseProduct(lambda);
	}

	// singleton + u no sampling:
	if(conditionalProbType == Conditions::STEM_TWO_EXT_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > dudt_hat = dxdt.getUnobservedNoSamplingStateProb();
		const Eigen::Ref< const Eigen::VectorXd > u_hat = x.getUnobservedNoSamplingStateProb();
		// duhatdt: step 1 - correction for phi since duhatdt -= (mu + lambda)*u_hat
		// duhatdt: step 1 - but we did duhatdt -= (phi + mu + lambda + delta)*u_hat
		// add mu too
		dudt_hat += mu + (phi + delta).cwiseProduct(u_hat);

		// dudt: step 3
		Eigen::MatrixXd resContractionVector(dudt_hat.rows(), 1);
		resContractionVector = computeTensorContractionVectorOpti<withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< const Eigen::VectorXd >, Eigen::Ref< const Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), u_hat, u_hat, t);
		dudt_hat += resContractionVector.cwiseProduct(lambda);

		// dsdt: step 1 - correction for phi since dsdt -= (mu + lambda + delta)*s
		// dsdt: step 1 - but we did dsdt -= (phi + mu + lambda + delta)*s
		Eigen::Ref< Eigen::VectorXd > dsdt_hat = dxdt.getSingletonNoSamplingStateProb();
		const Eigen::Ref< const Eigen::VectorXd > s_hat = x.getSingletonNoSamplingStateProb();

		// dsdt: step 1 - correction
		dsdt_hat += (phi + delta).cwiseProduct(s_hat);

		// dsdt: step 2
		resContractionVector = computeTensorContractionVectorOpti<withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< const Eigen::VectorXd >, Eigen::Ref< const Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), u_hat, s_hat, t);
		dsdt_hat += 2.*resContractionVector.cwiseProduct(lambda);

	}

}

} /* namespace Optimized */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */


#endif // LIKELIHOOD_KERNELS_CPU_OPTIMIZED_EIGENINTEGRATIONKERNEL_DEF_H_
