/*
 * EigenIntegrationKernel.cpp
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#include "EigenIntegrationKernel.h"

#include "../Misc/StiffnessRatio.h"
#ifndef LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENINTEGRATIONKERNEL_DEF_H_
#define LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENINTEGRATIONKERNEL_DEF_H_

#if defined(_OPENMP)

#include "EigenUtils.h"
#include "../EigenUtilsMatrix.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Likelihood {
namespace Kernels {
namespace CPU {
namespace OpenMP {

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::EigenIntegrationKernel(
			const size_t N_MAX_STATE_VECTOR,
			Phylogeny::Data::ContainerSharedPtr aPtrData,
			SynchronousEvents::ContainerSharedPtr aPtrSynchEventContainer,
			Tensor::ContainerSharedPtr aPtrTensorsContainer) :
					ptrData(aPtrData),
					ptrSynchEventContainer(aPtrSynchEventContainer),
					ptrTensorsContainer(aPtrTensorsContainer) {

	//assert(ptrTensorsContainer->getNumberOfState() > 2*Utils::Parallel::Manager::getInstance()->getMaxNThread());
	resFirstContractionU.resize(ptrTensorsContainer->getNumberOfState(), ptrTensorsContainer->getNumberOfState());

	isPrecomputedEtaAvailable = ptrTensorsContainer->getEta()->isContantThroughTime() &&
								ptrTensorsContainer->getMu()->isContantThroughTime() &&
								ptrTensorsContainer->getLambda()->isContantThroughTime() &&
								ptrTensorsContainer->getDelta()->isContantThroughTime() &&
								ptrTensorsContainer->getPhi()->isContantThroughTime() &&
								etaStructure != Tensor::ETA_QUASSE;

	if(isPrecomputedEtaAvailable) {
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
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::operator() ( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t ) {

	if(etaStructure == Tensor::ETA_DENSE && withCladoEvents && Utils::Parallel::Manager::getInstance()->useOpenMP()) {
		//doIntegrationStepChunk2GEMM(x, dxdt, t);
		doIntegrationStepBlock2GEMM(x, dxdt, t);
	} else {
		//doIntegrationStepChunk(x, dxdt, t);
		doIntegrationStepBlock(x, dxdt, t);
	}

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
std::pair<double, bool> EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::getStiffnessRatio(double t, const Eigen::VectorXd &u) const {
	StiffnessRatio sr(ptrTensorsContainer);
	return sr.estimateStiffnessRatio(t, u);
}

/*
template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::doIntegrationStepChunk( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t ) {


	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

	// Recover initial vectors
	const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);


	Eigen::MatrixXd vecSum;
	if(!isPrecomputedEtaAvailable){
		vecSum = phi + delta + mu + lambda;
	}


	const Eigen::Ref< const Eigen::VectorXd > u = x.getUnobservedStateProb();

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1")
	// For all probability vector: - (mu+delta+mu+lambda) * [u/p/...]+ eta*[u/p/...]
	{
		Eigen::Ref< Eigen::MatrixXd > dAlldt = dxdt.getStateProb();
		const Eigen::Ref< const Eigen::MatrixXd > all = x.getStateProb();

		size_t nCols = all.cols();
		if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {

			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
			std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
			size_t nThreads = chunks.size()-1;

			//omp_set_num_threads(nThreads);
			#pragma omp parallel
			{
				#pragma omp for schedule(static)
				for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
					Eigen::Ref< Eigen::MatrixXd > localDAlldt = dAlldt.block(0, chunks[iT], dAlldt.rows(), chunks[iT+1]-chunks[iT]);
					const Eigen::Ref< const Eigen::MatrixXd > localAll = all.block(0, chunks[iT], all.rows(), chunks[iT+1]-chunks[iT]);
					if(etaStructure == Tensor::ETA_SPARSE) {
						if(isPrecomputedEtaAvailable) {
							localDAlldt = precomputedEtaSparse*localAll;
						} else {
							const Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
							localDAlldt = eta[0]*localAll - (localAll.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*localAll;
						}
					} else if(etaStructure == Tensor::ETA_DENSE) {
						if(isPrecomputedEtaAvailable) {
							localDAlldt = precomputedEta*localAll;
						} else {
							const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
							localDAlldt = eta*localAll  - (localAll.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*localAll;
						}
					} else if(etaStructure == Tensor::ETA_QUASSE) {
						const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
						Eigen::MatrixXd vecSum = phi + delta + mu + lambda;
						localDAlldt =  defaultComputeQuasseEtaMatrix(eta(0,1), localAll) - (localAll.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
					}
				}
			}
		} else {
			if(etaStructure == Tensor::ETA_SPARSE) {
				if(isPrecomputedEtaAvailable) {
					dAlldt = precomputedEtaSparse*all;
				} else {
					const Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
					dAlldt = eta[0]*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
				}
			} else if(etaStructure == Tensor::ETA_DENSE) {
				if(isPrecomputedEtaAvailable) {
					dAlldt = precomputedEta*all;
				} else {
					const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
					dAlldt = eta*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
				}
			} else if(etaStructure == Tensor::ETA_QUASSE) {
				const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
				double scaledSigma = eta(0,1);
				dAlldt = defaultComputeQuasseEtaMatrix(scaledSigma, all) - (all.array().colwise()* vecSum.col(0).array()).matrix();
			}
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1")

	// For all [u/p]: d?dt += cst * lambda * omega * [p/u] u
	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")
	{
		Eigen::Ref< Eigen::MatrixXd > dupdt = dxdt.getUnobservedAndObservedStateProb();
		const Eigen::Ref< const Eigen::MatrixXd > up = x.getUnobservedAndObservedStateProb();

		if(withCladoEvents) {
			const Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
			computeFirstStepTensorContraction(omega, u, resFirstContractionU);
		} else {
			resFirstContractionU = u;
		}
		_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")

		_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_3")

		Eigen::MatrixXd resContractionMatrix(up.rows(), up.cols());
		Eigen::ArrayXd tmpArray1 = lambda.col(0).array();
		Eigen::ArrayXd tmpArray2 = 2.*Eigen::ArrayXd::Ones(dupdt.cols());
		tmpArray2(0) = 1.;

		size_t nCols = dupdt.cols();
		if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {

			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
			std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
			size_t nThreads = chunks.size()-1;

			//omp_set_num_threads(nThreads);
			#pragma omp parallel
			{
				#pragma omp for schedule(static)
				for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
					Eigen::Ref< Eigen::MatrixXd > localDupDt = dupdt.block(0, chunks[iT], dupdt.rows(), chunks[iT+1]-chunks[iT]);
					Eigen::Ref< Eigen::MatrixXd > localRCM = resContractionMatrix.block(0, chunks[iT], resContractionMatrix.rows(), chunks[iT+1]-chunks[iT]);
					const Eigen::Ref< const Eigen::ArrayXd > localArr2 = tmpArray2.segment(chunks[iT], chunks[iT+1]-chunks[iT]);

					if(!withCladoEvents) {
						localRCM = (up.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]).array().colwise() * resFirstContractionU.col(0).array()).matrix();
					} else {
						localRCM = resFirstContractionU * up.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]);
					}
					localDupDt += ((localRCM.array().colwise() * tmpArray1).rowwise() * localArr2.transpose()).matrix();
				}
			}
		} else {
			if(!withCladoEvents) {
				resContractionMatrix = (up.array().colwise() * resFirstContractionU.col(0).array()).matrix();
			} else {
				resContractionMatrix = resFirstContractionU * up;
			}
			dupdt += ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_3")

	// Unobserved: += mu
	{
		Eigen::Ref< Eigen::VectorXd > dudt = dxdt.getUnobservedStateProb();
		dudt += mu;
	}

	computeConditioningStateProb(x, dxdt, t);

	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

}
*/

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::doIntegrationStepBlock( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t ) {


	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

	// Recover initial vectors
	const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	const Eigen::Ref< const Eigen::VectorXd > u = x.getUnobservedStateProb();

	// Precompute
	Eigen::MatrixXd vecSum;
	if(!isPrecomputedEtaAvailable) {
		vecSum = phi + delta + mu + lambda;
	}


	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1")
	// For all probability vector: - (mu+delta+mu+lambda) * [u/p/...]+ eta*[u/p/...]
	{
		Eigen::Ref< Eigen::MatrixXd > dAlldt = dxdt.getStateProb();
		const Eigen::Ref< const Eigen::MatrixXd > all = x.getStateProb();

		if(etaStructure == Tensor::ETA_SPARSE) { // We don't know how to deal with block + sparse eta matrix
			size_t nCols = all.cols();
			if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {

				size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getNThread();
				std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
				size_t nThreads = chunks.size()-1;

				#pragma omp parallel
				{
					#pragma omp for schedule(static)
					for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
						Eigen::Ref< Eigen::MatrixXd > localDAlldt = dAlldt.block(0, chunks[iT], dAlldt.rows(), chunks[iT+1]-chunks[iT]);
						const Eigen::Ref< const Eigen::MatrixXd > localAll = all.block(0, chunks[iT], all.rows(), chunks[iT+1]-chunks[iT]);
						if(isPrecomputedEtaAvailable) {
							localDAlldt = precomputedEtaSparse*localAll;
						} else {
							const Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
							localDAlldt = eta[0]*localAll - (localAll.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*localAll;
						}
					}
				}
			} else {
				if(isPrecomputedEtaAvailable) {
					dAlldt = precomputedEtaSparse*all;
				} else {
					Tensor::sparseTensor_t &eta = ptrTensorsContainer->getEta()->getSparseTensor(t);
					dAlldt = eta[0]*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
				}
			}
		} else if(etaStructure == Tensor::ETA_DENSE) { // Eta is dense, block version
			size_t nRows = all.rows();
			size_t nCols = all.cols();
			_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "DEFINING_BLOCK_COST")
			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
			std::vector<Utils::Parallel::block_t> blocks(Utils::Parallel::Manager::getInstance()->defineOptimalBlocks(nAvailableThreads, nRows, nCols));
			size_t nThreads = blocks.size();
			_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "DEFINING_BLOCK_COST")

			if(nThreads > 1) {
				#pragma omp parallel
				{
					#pragma omp for schedule(static)
					for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
						Eigen::Ref< Eigen::MatrixXd > localDAlldt = dAlldt.block(blocks[iT].i, blocks[iT].j, blocks[iT].nI, blocks[iT].nJ);
						const Eigen::Ref< const Eigen::MatrixXd > localAll = all.block(blocks[iT].i, blocks[iT].j, blocks[iT].nI, blocks[iT].nJ);

						if(isPrecomputedEtaAvailable) {
								localDAlldt = precomputedEta.block(blocks[iT].i, 0, blocks[iT].nI, precomputedEta.cols())*all.block(0, blocks[iT].j, all.rows(), blocks[iT].nJ);
						} else {
								Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
								localDAlldt = eta.block(blocks[iT].i, 0, blocks[iT].nI, eta.cols())*all.block(0, blocks[iT].j, all.rows(), blocks[iT].nJ)
											  - (localAll.array().colwise() * vecSum.col(0).segment(blocks[iT].i, blocks[iT].nI).array()).matrix();
												//vecSum.col(0).segment(blocks[iT].i, blocks[iT].nI).asDiagonal() * localAll;
						}
					}
				}
			} else {
				if(isPrecomputedEtaAvailable) {
					dAlldt = precomputedEta*all;
				} else {
					Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
					dAlldt = eta*all - (all.array().colwise() * vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
				}
			}
		} else if(etaStructure == Tensor::ETA_QUASSE) { // No block here
			size_t nCols = all.cols();
			if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {

				size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
				std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
				size_t nThreads = chunks.size()-1;

				#pragma omp parallel
				{
					#pragma omp for schedule(static)
					for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
						Eigen::Ref< Eigen::MatrixXd > localDAlldt = dAlldt.block(0, chunks[iT], dAlldt.rows(), chunks[iT+1]-chunks[iT]);
						const Eigen::Ref< const Eigen::MatrixXd > localAll = all.block(0, chunks[iT], all.rows(), chunks[iT+1]-chunks[iT]);

						const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
						double scaledSigma = eta(0,1);
						localDAlldt = defaultComputeQuasseEtaMatrix(scaledSigma, localAll) - (localAll.array().colwise()* vecSum.col(0).array()).matrix();

					}
				}
			} else {
				const Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
				double scaledSigma = eta(0,1);
				dAlldt =  defaultComputeQuasseEtaMatrix(scaledSigma, all) - (all.array().colwise()* vecSum.col(0).array()).matrix();
			}
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1")

	// For all [u/p]: d?dt += cst * lambda * omega * [p/u] u
	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")
	{
		Eigen::Ref< Eigen::MatrixXd > dupdt = dxdt.getUnobservedAndObservedStateProb();
		const Eigen::Ref< const Eigen::MatrixXd > up = x.getUnobservedAndObservedStateProb();

		if(withCladoEvents) {
			const Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
			computeFirstStepTensorContraction(omega, u, resFirstContractionU);
		} else {
			resFirstContractionU = u;
		}
		_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")

		_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_3")

		Eigen::MatrixXd resContractionMatrix(up.rows(), up.cols());
		Eigen::ArrayXd tmpArray1 = lambda.col(0).array();
		Eigen::ArrayXd tmpArray2 = 2.*Eigen::ArrayXd::Ones(dupdt.cols());
		tmpArray2(0) = 1.;

		if(!withCladoEvents) { // column -- we need column chunks
			size_t nCols = dupdt.cols();
			if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {

				size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
				std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
				size_t nThreads = chunks.size()-1;

				//omp_set_num_threads(nThreads);
				#pragma omp parallel
				{
					#pragma omp for schedule(static)
					for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
						Eigen::Ref< Eigen::MatrixXd > localDupDt = dupdt.block(0, chunks[iT], dupdt.rows(), chunks[iT+1]-chunks[iT]);
						Eigen::Ref< Eigen::MatrixXd > localRCM = resContractionMatrix.block(0, chunks[iT], resContractionMatrix.rows(), chunks[iT+1]-chunks[iT]);
						const Eigen::Ref< const Eigen::ArrayXd > localArr2 = tmpArray2.segment(chunks[iT], chunks[iT+1]-chunks[iT]);

						localRCM = (up.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]).array().colwise() * resFirstContractionU.col(0).array()).matrix();
						localDupDt += ((localRCM.array().colwise() * tmpArray1).rowwise() * localArr2.transpose()).matrix();
					}
				}
			} else {
				resContractionMatrix = (up.array().colwise() * resFirstContractionU.col(0).array()).matrix();
				dupdt += ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
			}
		} else { // Matrix (first contraction phase) we need blocks
			size_t nRows = dupdt.rows();
			size_t nCols = dupdt.cols();
			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
			std::vector<Utils::Parallel::block_t> blocks(Utils::Parallel::Manager::getInstance()->defineOptimalBlocks(nAvailableThreads, nRows, nCols));
			size_t nThreads = blocks.size();

			if(nThreads > 1) {
				//omp_set_num_threads(nThreads);
				#pragma omp parallel
				{
					#pragma omp for schedule(static)
					for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
						Eigen::Ref< Eigen::MatrixXd > localDupDt = dupdt.block(blocks[iT].i, blocks[iT].j, blocks[iT].nI, blocks[iT].nJ);
						Eigen::Ref< Eigen::MatrixXd > localRCM = resContractionMatrix.block(blocks[iT].i, blocks[iT].j, blocks[iT].nI, blocks[iT].nJ);
						const Eigen::Ref< const Eigen::ArrayXd > localArr1 = tmpArray1.segment(blocks[iT].i, blocks[iT].nI);
						const Eigen::Ref< const Eigen::ArrayXd > localArr2 = tmpArray2.segment(blocks[iT].j, blocks[iT].nJ);

						localRCM = resFirstContractionU.block(blocks[iT].i, 0, blocks[iT].nI, resFirstContractionU.cols()) * up.block(0, blocks[iT].j, up.rows(), blocks[iT].nJ);
						localDupDt += ((localRCM.array().colwise() * localArr1).rowwise() * localArr2.transpose()).matrix();
					}
				}
			} else {
				//std::cout << "Matrix resContr at t = " <<  t << " : " << std::endl << resFirstContractionU << std::endl << "--------------------------------------" << std::endl;
				resContractionMatrix = resFirstContractionU * up;
				dupdt += ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
			}
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_3")



	// Unobserved: += mu
	{
		Eigen::Ref< Eigen::VectorXd > dudt = dxdt.getUnobservedStateProb();
		dudt += mu;
	}

	computeConditioningStateProb(x, dxdt, t);

	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

}

/*
template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::doIntegrationStepChunk2GEMM( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t ) {


	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

	// Recover initial vectors
	const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	const Eigen::Ref< const Eigen::VectorXd > u = x.getUnobservedStateProb();

	// Precompute
	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")

	if(withCladoEvents) {
		const Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		computeFirstStepTensorContraction(omega, u, resFirstContractionU);
	} else {
		resFirstContractionU = u;
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")


	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1+3")

	size_t nThreadPhase1 = std::ceil(Utils::Parallel::Manager::getInstance()->getMaxNThread()/2.);
	size_t nThreadPhase2 = std::floor(Utils::Parallel::Manager::getInstance()->getMaxNThread()/2.);
	//std::cout << "nThreadPhase 1 = " << nThreadPhase1 << " -- nThreadPhase 2 = " << nThreadPhase2 << std::endl;

	// REQUIRED VARIABLES
	Eigen::MatrixXd vecSum;
	if(!isPrecomputedEtaAvailable) {
		vecSum = phi + delta + mu + lambda;
	}

	Eigen::Ref< Eigen::MatrixXd > dAlldt = dxdt.getStateProb();
	const Eigen::Ref< const Eigen::MatrixXd > all = x.getStateProb();
	const Eigen::Ref< const Eigen::MatrixXd > up = x.getUnobservedAndObservedStateProb();

	Eigen::MatrixXd resContractionMatrix(up.rows(), up.cols());
	Eigen::ArrayXd tmpArray1 = lambda.col(0).array();
	Eigen::ArrayXd tmpArray2 = 2.*Eigen::ArrayXd::Ones(up.cols());
	tmpArray2(0) = 1.;


	// GEMM 1 decomposition
	size_t nColsGEMM1 = all.cols();
	bool parallelGEMM1 = Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nColsGEMM1, nThreadPhase1);
	size_t nThreadGEMM1 = 1;
	std::vector<size_t> chunksGEMM1;
	if(parallelGEMM1) {
		chunksGEMM1 = Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nThreadPhase1, nColsGEMM1);
		nThreadGEMM1 = chunksGEMM1.size()-1;
	}
	//std::cout << "GEMM1 nThread = " << nThreadGEMM1 << " :: bCols = " << nColsGEMM1 << std::endl;


	// GEMM 2 decomposition
	size_t nColsGEMM2 = up.cols();
	bool parallelGEMM2 = Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nColsGEMM2, nThreadPhase2);
	size_t nThreadGEMM2 = 1;
	std::vector<size_t> chunksGEMM2;
	if(parallelGEMM2) {
		chunksGEMM2 = Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nThreadPhase2, nColsGEMM2);
		nThreadGEMM2 = chunksGEMM2.size()-1;
	}
	//std::cout << "GEMM2 nThread = " << nThreadGEMM2 << " :: bCols = " << nColsGEMM2 << std::endl;

	// Total threads
	size_t nTotalThreads = nThreadGEMM1+nThreadGEMM2;

	//omp_set_num_threads(nTotalThreads);
	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for(size_t iT=0; iT<nTotalThreads; ++iT) { // Potential for pragma openmp

			if(iT < nThreadGEMM1) {
				if(parallelGEMM1) {
					//std::cout << "GEMM1 " << iT << " -- chunk : " << chunksGEMM1[iT] << ", " << chunksGEMM1[iT+1] << ", nCols = " << dAlldt.cols() << std::endl;
					Eigen::Ref< Eigen::MatrixXd > localDAlldt = dAlldt.block(0, chunksGEMM1[iT], dAlldt.rows(), chunksGEMM1[iT+1]-chunksGEMM1[iT]);
					const Eigen::Ref< const Eigen::MatrixXd > localAll = all.block(0, chunksGEMM1[iT], all.rows(), chunksGEMM1[iT+1]-chunksGEMM1[iT]);

					if(isPrecomputedEtaAvailable) {
						localDAlldt = precomputedEta*localAll;
					} else {
						Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
						localDAlldt = eta*localAll - (localAll.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*localAll;
					}

				} else {
					//std::cout << "GEMM1 " << iT << " SEQ " << std::endl;
					assert(nThreadGEMM1 == 1);

					if(isPrecomputedEtaAvailable) {
						dAlldt = precomputedEta*all;
					} else {
						Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
						dAlldt = eta*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
					}
				}

			} else {

				size_t iT2 = iT-nThreadGEMM1;

				if(parallelGEMM2) {
					//std::cout << "GEMM2 " << iT << " :: " << iT2 << " -- chunk : " << chunksGEMM2[iT2] << ", " << chunksGEMM2[iT2+1] << ", nCols = " << up.cols() << std::endl;
					//std::cout << resContractionMatrix.rows() << " x " << resContractionMatrix.cols() << " :: " << up.rows() << " x " << up.cols() << std::endl;
					Eigen::Ref< Eigen::MatrixXd > localRCM = resContractionMatrix.block(0, chunksGEMM2[iT2], up.rows(), chunksGEMM2[iT2+1]-chunksGEMM2[iT2]);
					const Eigen::Ref< const Eigen::ArrayXd > localArr2 = tmpArray2.segment(chunksGEMM2[iT2], chunksGEMM2[iT2+1]-chunksGEMM2[iT2]);

					localRCM = resFirstContractionU * up.block(0, chunksGEMM2[iT2], up.rows(), chunksGEMM2[iT2+1]-chunksGEMM2[iT2]);
					localRCM = ((localRCM.array().colwise() * tmpArray1).rowwise() * localArr2.transpose()).matrix();

				} else {
					//std::cout << "GEMM2 " << iT << " SEQ " << std::endl;
					assert(nThreadGEMM2 == 1);
					resContractionMatrix = resFirstContractionU * up;
					resContractionMatrix = ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
				}
			}
		}
	}

	{
		size_t nCols = all.cols();
		Eigen::Ref< Eigen::MatrixXd > dupdt = dxdt.getUnobservedAndObservedStateProb();

		if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {
			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
			std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
			size_t nThreads = chunks.size()-1;

			//omp_set_num_threads(nThreads);
			#pragma omp parallel
			{
				#pragma omp for schedule(static)
				for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
					dupdt.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]) += resContractionMatrix.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]);
				}
			}
		} else {
			dupdt += resContractionMatrix;
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1+3")


	// Unobserved: += mu
	{
		Eigen::Ref< Eigen::VectorXd > dudt = dxdt.getUnobservedStateProb();
		dudt += mu;
	}

	computeConditioningStateProb(x, dxdt, t);

	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

}
*/

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::doIntegrationStepBlock2GEMM( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t ) {

	assert(etaStructure == Tensor::ETA_DENSE);

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

	// Recover initial vectors
	const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	const Eigen::Ref< const Eigen::VectorXd > u = x.getUnobservedStateProb();

	// Precompute
	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")

	if(withCladoEvents) {
		const Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
		computeFirstStepTensorContraction(omega, u, resFirstContractionU);
	} else {
		resFirstContractionU = u;
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_2")


	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1+3")

	size_t nThreadPhase1 = std::ceil(Utils::Parallel::Manager::getInstance()->getNThread()/2.);
	size_t nThreadPhase2 = std::floor(Utils::Parallel::Manager::getInstance()->getNThread()/2.);
	//std::cout << "nThreadPhase 1 = " << nThreadPhase1 << " -- nThreadPhase 2 = " << nThreadPhase2 << std::endl;

	// REQUIRED VARIABLES
	Eigen::MatrixXd vecSum = phi + delta + mu + lambda;
	Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);

	Eigen::Ref< Eigen::MatrixXd > dAlldt = dxdt.getStateProb();
	const Eigen::Ref< const Eigen::MatrixXd > all = x.getStateProb();
	const Eigen::Ref< const Eigen::MatrixXd > up = x.getUnobservedAndObservedStateProb();

	Eigen::MatrixXd resContractionMatrix(up.rows(), up.cols());
	Eigen::ArrayXd tmpArray1 = lambda.col(0).array();
	Eigen::ArrayXd tmpArray2 = 2.*Eigen::ArrayXd::Ones(up.cols());
	tmpArray2(0) = 1.;


	// GEMM 1 decomposition
	size_t nRowsGEMM1 = all.rows();
	size_t nColsGEMM1 = all.cols();
	size_t nAvailableThreadsGEMM1 = nThreadPhase1;
	std::vector<Utils::Parallel::block_t> blocksGEMM1(Utils::Parallel::Manager::getInstance()->defineOptimalBlocks(nAvailableThreadsGEMM1, nRowsGEMM1, nColsGEMM1));
	size_t nThreadGEMM1 = blocksGEMM1.size();
	//std::cout << "GEMM1 nThread = " << nThreadGEMM1 << " :: bCols = " << nColsGEMM1 << std::endl;


	// GEMM 2 decomposition
	size_t nRowsGEMM2 = up.rows();
	size_t nColsGEMM2 = up.cols();
	size_t nAvailableThreadsGEMM2 = nThreadPhase2;
	std::vector<Utils::Parallel::block_t> blocksGEMM2(Utils::Parallel::Manager::getInstance()->defineOptimalBlocks(nAvailableThreadsGEMM2, nRowsGEMM2, nColsGEMM2));
	size_t nThreadGEMM2 = blocksGEMM2.size();

	//std::cout << "GEMM2 nThread = " << nThreadGEMM2 << " :: bCols = " << nColsGEMM2 << std::endl;

	// Total threads
	size_t nTotalThreads = nThreadGEMM1+nThreadGEMM2;

	//omp_set_num_threads(nTotalThreads);
	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for(size_t iT=0; iT<nTotalThreads; ++iT) { // Potential for pragma openmp

			if(iT < nThreadGEMM1) {
				if(nThreadGEMM1 > 1) {
					Eigen::Ref< Eigen::MatrixXd > localDAlldt = dAlldt.block(blocksGEMM1[iT].i, blocksGEMM1[iT].j, blocksGEMM1[iT].nI, blocksGEMM1[iT].nJ);
					const Eigen::Ref< const Eigen::MatrixXd > localAll = all.block(blocksGEMM1[iT].i, blocksGEMM1[iT].j, blocksGEMM1[iT].nI, blocksGEMM1[iT].nJ);
					if(isPrecomputedEtaAvailable) {
						localDAlldt =  precomputedEta.block(blocksGEMM1[iT].i, 0, blocksGEMM1[iT].nI, eta.cols())*all.block(0, blocksGEMM1[iT].j, all.rows(), blocksGEMM1[iT].nJ);
					} else {
						localDAlldt =  eta.block(blocksGEMM1[iT].i, 0, blocksGEMM1[iT].nI, eta.cols())*all.block(0, blocksGEMM1[iT].j, all.rows(), blocksGEMM1[iT].nJ)
									  - (localAll.array().colwise() * vecSum.col(0).segment(blocksGEMM1[iT].i, blocksGEMM1[iT].nI).array()).matrix();
										//vecSum.col(0).segment(blocksGEMM1[iT].i, blocksGEMM1[iT].nI).asDiagonal() * localAll;
					}
				} else {
					if(isPrecomputedEtaAvailable) {
						dAlldt = precomputedEta*all;
					} else {
						Eigen::MatrixXd &eta = ptrTensorsContainer->getEigenMatrixEta(t);
						dAlldt = eta*all - (all.array().colwise()* vecSum.col(0).array()).matrix();//vecSum.col(0).asDiagonal()*all;
					}
				}
			} else {

				size_t iT2 = iT-nThreadGEMM2;

				if(nThreadGEMM2 > 1) {
					Eigen::Ref< Eigen::MatrixXd > localRCM = resContractionMatrix.block(blocksGEMM2[iT2].i, blocksGEMM2[iT2].j, blocksGEMM2[iT2].nI, blocksGEMM2[iT2].nJ);
					const Eigen::Ref< const Eigen::ArrayXd > localArr1 = tmpArray1.segment(blocksGEMM2[iT2].i, blocksGEMM2[iT2].nI);
					const Eigen::Ref< const Eigen::ArrayXd > localArr2 = tmpArray2.segment(blocksGEMM2[iT2].j, blocksGEMM2[iT2].nJ);

					localRCM = resFirstContractionU.block(blocksGEMM2[iT2].i, 0, blocksGEMM2[iT2].nI, resFirstContractionU.cols()) * up.block(0, blocksGEMM2[iT2].j, up.rows(), blocksGEMM2[iT2].nJ);
					localRCM = ((localRCM.array().colwise() * localArr1).rowwise() * localArr2.transpose()).matrix();
				} else {
					assert(nThreadGEMM2 == 1);
					resContractionMatrix = resFirstContractionU * up;
					resContractionMatrix = ((resContractionMatrix.array().colwise() * tmpArray1).rowwise() * tmpArray2.transpose()).matrix();
				}
			}
		}
	}

	{
		Eigen::Ref< Eigen::MatrixXd > dupdt = dxdt.getUnobservedAndObservedStateProb();
		size_t nCols = dupdt.cols();

		if(Utils::Parallel::Manager::getInstance()->useBlockParallelOMP(nCols)) {
			size_t nAvailableThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
			std::vector<size_t> chunks(Utils::Parallel::Manager::getInstance()->defineOptimalChunks(nAvailableThreads, nCols));
			size_t nThreads = chunks.size()-1;

			/*std::cout << "nThreads = " << nThreads << std::endl;
			for(size_t iC=0; iC<chunks.size(); ++iC) {
				std::cout << iC << " :: chunks = " << chunks[iC] << std::endl;
			}
			std::cout << dupdt.rows() << " -- " << dupdt.cols() << "  :: " << resContractionMatrix.rows() << " -- " << resContractionMatrix.cols() << std::endl;
			 */

			//omp_set_num_threads(nThreads);
			#pragma omp parallel
			{
				#pragma omp for schedule(static)
				for(size_t iT=0; iT<nThreads; ++iT) { // Potential for pragma openmp
					//std::cout << "iT = " << iT << std::endl;
					//std::cout << "Chunks : " << chunks[iT] << " -- " << chunks[iT+1] << std::endl;
					dupdt.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]) += resContractionMatrix.block(0, chunks[iT], up.rows(), chunks[iT+1]-chunks[iT]);
				}
			}
		} else {
			dupdt += resContractionMatrix;
		}
	}
	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP_1+3")


	// Unobserved: += mu
	{
		Eigen::Ref< Eigen::VectorXd > dudt = dxdt.getUnobservedStateProb();
		dudt += mu;
	}

	computeConditioningStateProb(x, dxdt, t);

	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INTEGRATION_STEP")

}

template <  Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void EigenIntegrationKernel<conditionalProbType, etaStructure, withCladoEvents>::computeConditioningStateProb( const Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &x , Likelihood::StateType::OpenMP::EigenState<conditionalProbType> &dxdt , double t ) {

	const Eigen::MatrixXd &mu = ptrTensorsContainer->getEigenVecMu(t);
	const Eigen::MatrixXd &lambda = ptrTensorsContainer->getEigenVecLambda(t);
	const Eigen::MatrixXd &phi = ptrTensorsContainer->getEigenVecPhi(t);
	const Eigen::MatrixXd &delta = ptrTensorsContainer->getEigenVecDelta(t);

	const Eigen::Ref< const Eigen::VectorXd > u = x.getUnobservedStateProb();

	// singleton:
	if(conditionalProbType == Conditions::STEM_TWO_SAMPLES) {
		Eigen::Ref< Eigen::VectorXd > dsdt = dxdt.getSingletonStateProb();
		const Eigen::Ref< const Eigen::VectorXd > s = x.getSingletonStateProb();
		// dsdt: step 1 - correction for phi since dsdt -= (mu + lambda + delta)*s
		// dsdt: step 1 - but we did dsdt -= (phi + mu + lambda + delta)*s
//		dsdt += phi.cwiseProduct(s);
		// dsdt: step 2
		Eigen::MatrixXd resContractionVector;
		computeSecondStepTensorContractionVector< withCladoEvents,Eigen::Ref< const Eigen::VectorXd > >(resFirstContractionU, s, resContractionVector);
		dsdt += 2.0 * resContractionVector.cwiseProduct(lambda);

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
		Eigen::VectorXd resContraction = computeTensorContractionVectorOpti< withCladoEvents, Tensor::BaseTensorSharedPtr, Eigen::Ref< const Eigen::VectorXd >, Eigen::Ref< const Eigen::VectorXd > >(ptrTensorsContainer->getOmega(), u_hat, u_hat, t);
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
		if(withCladoEvents) {
			Tensor::sparseTensor_t &omega = ptrTensorsContainer->getOmega()->getSparseTensor(t);
			computeFirstStepTensorContraction(omega, u_hat, resFirstContractionU);
		} else {
			resFirstContractionU = u_hat;
		}
		computeSecondStepTensorContractionVector< withCladoEvents,Eigen::Ref< const Eigen::VectorXd > >(resFirstContractionU, u_hat, resContractionVector);
		dudt_hat += resContractionVector.cwiseProduct(lambda);

		// dsdt: step 1 - correction for phi since dsdt -= (mu + lambda + delta)*s
		// dsdt: step 1 - but we did dsdt -= (phi + mu + lambda + delta)*s
		Eigen::Ref< Eigen::VectorXd > dsdt_hat = dxdt.getSingletonNoSamplingStateProb();
		const Eigen::Ref< const Eigen::VectorXd > s_hat = x.getSingletonNoSamplingStateProb();

		// dsdt: step 1 - correction
		dsdt_hat += (phi + delta).cwiseProduct(s_hat);

		// dsdt: step 2
		computeSecondStepTensorContractionVector< withCladoEvents,Eigen::Ref< const Eigen::VectorXd > >(resFirstContractionU, s_hat, resContractionVector);
		dsdt_hat += 2.*resContractionVector.cwiseProduct(lambda);
	}
}


} /* namespace OpenMP */
} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif //defined(_OPENMP)
#endif // LIKELIHOOD_KERNELS_CPU_OPENMP_EIGENINTEGRATIONKERNEL_DEF_H_
