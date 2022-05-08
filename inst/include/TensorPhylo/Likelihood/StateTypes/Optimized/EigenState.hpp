/*
 * EigenState.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATE_HPP_
#define LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATE_HPP_

#include "EigenState.h"

#include "Utils/MemoryPool/EigenCPU.h"
#include "../Utils.h"

namespace Likelihood {
namespace StateType {
namespace Optimized {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::idSeq = 0;

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>::EigenState() : id(idSeq++), N_VECTOR(Likelihood::Conditions::getNVector(conditionalProbType)) {
	stateProb = NULL;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>::EigenState(const size_t aNEdges) : id(idSeq++), N_VECTOR(Likelihood::Conditions::getNVector(conditionalProbType)) {
	stateProb = NULL;
	initMemory();
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>::EigenState(const EigenState &aEigenState) :
		id(idSeq++), N_VECTOR(Likelihood::Conditions::getNVector(conditionalProbType)),
		vecScaling(aEigenState.vecScaling) {

	if(aEigenState.stateProb != NULL) {
		initMemory();
		size_t nStates = aEigenState.stateProb->matrix.rows();
		size_t nCols = N_VECTOR + aEigenState.vecProbToEdgeMapping.size();
		stateProb->matrix.block(0, 0, nStates, nCols) = aEigenState.stateProb->matrix.block(0, 0, nStates, nCols);
		vecProbToEdgeMapping = aEigenState.vecProbToEdgeMapping;
		vecScaling = aEigenState.vecScaling;
	}  else {
		stateProb = NULL;
	}
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>::~EigenState() {
	//std::cout << " Delete ID = " << id << std::endl;

	if(stateProb != NULL) {
		Utils::MemoryPool::eigenCPU().freeMatrix(stateProb);
		stateProb = NULL;
	}

}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::setVecProbToEdgeMapping(std::vector<int> &aMapping) {
	vecProbToEdgeMapping = aMapping;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
std::vector<int>& EigenState<conditionalProbType>::getVecProbToEdgeMapping() {
	return vecProbToEdgeMapping;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const std::vector<int>& EigenState<conditionalProbType>::getVecProbToEdgeMapping() const {
	return vecProbToEdgeMapping;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::allocateVecProbForEdge(size_t idEdge) {

	if(stateProb == NULL) {
		initMemory();
	}

	assert(vecProbToEdgeMapping.size() < (size_t)(stateProb->matrix.cols()-N_VECTOR));

	int idVec = vecProbToEdgeMapping.size();

	vecProbToEdgeMapping.push_back(idEdge);
	vecScaling.push_back(0.);

	return idVec;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::removeVecProbForEdge(size_t idEdge) {
	// Look at the vec prob position
	std::vector<int>::iterator itFind = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), idEdge);
	assert(itFind != vecProbToEdgeMapping.end());

	// Get the position
	size_t iVecProb = std::distance(vecProbToEdgeMapping.begin(), itFind);

	Eigen::Ref< Eigen::MatrixXd > observedStateProb = getObservedStateProb();

	// If its not the last probability vector: swap
	size_t idLastCol = vecProbToEdgeMapping.size()-1;
	if(iVecProb != idLastCol) {
		std::swap(vecScaling[iVecProb], vecScaling[idLastCol]);
		std::swap(vecProbToEdgeMapping[iVecProb], vecProbToEdgeMapping[idLastCol]);

		observedStateProb.col(iVecProb).swap(observedStateProb.col(idLastCol));
		iVecProb = idLastCol;
	}

	// Set to zero in matrix
	observedStateProb.col(idLastCol).setZero();

	// Erase in mapping
	vecProbToEdgeMapping.erase(vecProbToEdgeMapping.begin()+idLastCol);
	vecScaling.erase(vecScaling.begin()+idLastCol);
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::updateVecProbForEdgeMapping(size_t fromIdEdge, size_t toIdEdge) {
	// Look at the vec prob position
	std::vector<int>::iterator itFind = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), fromIdEdge);
	assert(itFind != vecProbToEdgeMapping.end());

	// Get the position (could directly use the iterator but illisible)
	size_t iVecProb = std::distance(vecProbToEdgeMapping.begin(), itFind);

	vecProbToEdgeMapping[iVecProb]= toIdEdge;
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::findVecProbIdForEdgeId(size_t iE) const {
	// Look at the vec prob position
	std::vector<int>::const_iterator itFind = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), iE);
	assert(itFind != vecProbToEdgeMapping.end());

	// Get the position (could directly use the iterator but illisible)
	size_t iVecProb = std::distance(vecProbToEdgeMapping.begin(), itFind);
	return iVecProb;
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::resize(size_t aNEdges) {
	if(aNEdges == 0) {
		vecProbToEdgeMapping.clear();
		vecScaling.clear();
		return;
	}

	if(stateProb == NULL) {
		initMemory();
	}

	size_t initSize = size();
	if(initSize < aNEdges) {
		for(size_t iE=initSize; iE<aNEdges; ++iE) {
			vecProbToEdgeMapping.push_back(-1);
			vecScaling.push_back(0.);
		}
	} else if(initSize > aNEdges) {
		for(size_t iE=aNEdges; iE<initSize; ++iE) {
			vecProbToEdgeMapping.pop_back();
			vecScaling.pop_back();
		}
	}

	assert(size() == aNEdges);

}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::getNVector() const {
	return N_VECTOR;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::VectorXd > EigenState<conditionalProbType>::getUnobservedStateProb() {
	if(stateProb == NULL) {
		initMemory();
	}
	return stateProb->matrix.col(getIdxU());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::VectorXd > EigenState<conditionalProbType>::getUnobservedNoSamplingStateProb() {
	if(stateProb == NULL) {
		initMemory();
	}
	return stateProb->matrix.col(getIdxUHat());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::VectorXd > EigenState<conditionalProbType>::getSingletonStateProb() {
	if(stateProb == NULL) {
		initMemory();
	}
	return stateProb->matrix.col(getIdxS());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::VectorXd > EigenState<conditionalProbType>::getSingletonNoSamplingStateProb() {
	if(stateProb == NULL) {
		initMemory();
	}
	return stateProb->matrix.col(getIdxSHat());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::VectorXd > EigenState<conditionalProbType>::getObservedStateProb(size_t iVec) {
	if(stateProb == NULL) {
		initMemory();
	}
	return stateProb->matrix.col(N_VECTOR+iVec);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::MatrixXd > EigenState<conditionalProbType>::getObservedStateProb() {
	if(stateProb == NULL) {
		initMemory();
	}
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = vecProbToEdgeMapping.size();
	return stateProb->matrix.block(0, N_VECTOR, nStates, nCols);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::MatrixXd > EigenState<conditionalProbType>::getStateProb() {
	if(stateProb == NULL) {
		initMemory();
	}
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = N_VECTOR + vecProbToEdgeMapping.size();
	return stateProb->matrix.block(0, 0, nStates, nCols);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Eigen::Ref< Eigen::MatrixXd > EigenState<conditionalProbType>::getUnobservedAndObservedStateProb() {
	assert(stateProb != NULL);
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = 1 + vecProbToEdgeMapping.size();
	return stateProb->matrix.block(0, N_VECTOR-1, nStates, nCols);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::VectorXd > EigenState<conditionalProbType>::getUnobservedStateProb() const {
	assert(stateProb != NULL);
	return stateProb->matrix.col(getIdxU());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::VectorXd > EigenState<conditionalProbType>::getUnobservedNoSamplingStateProb() const {
	assert(stateProb != NULL);
	return stateProb->matrix.col(getIdxUHat());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::VectorXd > EigenState<conditionalProbType>::getSingletonStateProb() const {
	assert(stateProb != NULL);
	return stateProb->matrix.col(getIdxS());
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::VectorXd > EigenState<conditionalProbType>::getSingletonNoSamplingStateProb() const {
	assert(stateProb != NULL);
	return stateProb->matrix.col(getIdxSHat());
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::MatrixXd > EigenState<conditionalProbType>::getObservedStateProb() const {
	assert(stateProb != NULL);
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = vecProbToEdgeMapping.size();
	return stateProb->matrix.block(0, N_VECTOR, nStates, nCols);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::MatrixXd > EigenState<conditionalProbType>::getStateProb() const {
	assert(stateProb != NULL);
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = N_VECTOR + vecProbToEdgeMapping.size();
	return stateProb->matrix.block(0, 0, nStates, nCols);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
const Eigen::Ref< const Eigen::MatrixXd > EigenState<conditionalProbType>::getUnobservedAndObservedStateProb() const {
	assert(stateProb != NULL);
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = 1 + vecProbToEdgeMapping.size();
	return stateProb->matrix.block(0, N_VECTOR-1, nStates, nCols);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::size() const {
	if(stateProb == NULL) {
		return 0;
	} else {
		return vecProbToEdgeMapping.size();
	}
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
double EigenState<conditionalProbType>::defineNormInf() const {
	//std::cout << " Norm inf of id  = " << id << std::endl;
	assert(stateProb != NULL);
	const Eigen::Ref< const Eigen::MatrixXd > blockSProb = getStateProb();
    double absMax = blockSProb.lpNorm<Eigen::Infinity>();
    return absMax;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::addMult(double factor, const EigenState &state ) {
	assert(stateProb != NULL);
	assert(state.size() == size());

	getStateProb() += factor*state.getStateProb();
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::initMult(double factor, const EigenState &state ) {
	resize(state.size());

	assert(stateProb != NULL);
	getStateProb() = factor * state.getStateProb();

	vecProbToEdgeMapping = state.vecProbToEdgeMapping;
	vecScaling = state.vecScaling;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::odeIntRelativeError(double m_eps_abs, double m_eps_rel, double m_a_x, double m_a_dxdt, const EigenState &state1, const EigenState &state2) {
	assert(state1.size() == size() && state2.size() == size());

	assert(stateProb != NULL);
	Eigen::Ref< Eigen::MatrixXd > blockSProb = getStateProb();
	Eigen::MatrixXd denom = (m_eps_abs + (m_eps_rel * ( m_a_x * state1.getStateProb() + m_a_dxdt * state2.getStateProb().cwiseAbs())).array()).matrix();
	blockSProb = blockSProb.cwiseAbs().cwiseQuotient(denom);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>& EigenState<conditionalProbType>::operator+=( const double &val ) {
	//std::cout << " ADD val = " << val << " to id = " << id << std::endl;

	assert(stateProb != NULL);
	getStateProb() = (getStateProb().array()  + val).matrix();

	return *this;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>& EigenState<conditionalProbType>::operator+=( const EigenState &state ) {

	assert(state.size() == size());
	vecProbToEdgeMapping = state.vecProbToEdgeMapping;
	vecScaling = state.vecScaling;

	getStateProb() += state.getStateProb();

	return *this;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>& EigenState<conditionalProbType>::operator*=( const double a ) {
	//std::cout << " Multiply id = " << id << " by = " << a << std::endl;

	assert(stateProb != NULL);
	getStateProb() *= a;
	return *this;
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType>& EigenState<conditionalProbType>::operator=(const EigenState &state) {
	if(	state.stateProb != NULL) {
		resize(state.size());

		vecProbToEdgeMapping = state.vecProbToEdgeMapping;
		vecScaling = state.vecScaling;
		getStateProb() = state.getStateProb();

	} else  if(stateProb != NULL) {
			Utils::MemoryPool::eigenCPU().freeMatrix(stateProb);
			stateProb = NULL;
			vecProbToEdgeMapping.clear();
			vecScaling.clear();
	}

	return *this;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::initMemory() {
	stateProb = Utils::MemoryPool::eigenCPU().allocateMatrix();
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::getIdxU() const {
	return N_VECTOR-1;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::getIdxUHat() const {
	if(conditionalProbType == Likelihood::Conditions::ROOT_SURVIVAL ||
	   conditionalProbType == Likelihood::Conditions::ROOT_MRCA ||
	   conditionalProbType == Likelihood::Conditions::STEM_SURVIVAL ||
	   conditionalProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
		return 0;
	}
	assert(false);
	return -1;
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::getIdxS() const {
	if(conditionalProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
		return 0;
	}

	assert(false);
	return -1;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
size_t EigenState<conditionalProbType>::getIdxSHat() const {
	if(conditionalProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
		return 1;
	}
	assert(false);
	return -1;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::roundNegativeProbabilityToZero() {

	assert(stateProb != NULL);
	size_t nStates = stateProb->matrix.rows();
	size_t nCols = N_VECTOR + vecProbToEdgeMapping.size();
	if((stateProb->matrix.block(0, 0, nStates, nCols).array() < 0.).any()) {
		for(size_t iR=0; iR<nStates; ++iR) {
			for(size_t iC=0; iC<nCols; ++iC) {
				if(stateProb->matrix(iR, iC) < 0.) {
					stateProb->matrix(iR, iC) = 0.;
				}
			}
		}
	}

}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::rescaleAll() {

	if(stateProb == NULL) return;

	size_t nCols = vecProbToEdgeMapping.size();
	for(size_t iC=0; iC<nCols; ++iC) {
		vecScaling[iC] += rescaleProbabilityVector(getObservedStateProb(iC));
	}


}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
double EigenState<conditionalProbType>::getScalingFactorByVecPos(size_t iV) const {
	assert(iV < vecScaling.size());
	return vecScaling[iV];
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
double EigenState<conditionalProbType>::getScalingFactorByEdgeId(size_t iE) const {
	// Look at the vec prob position
	std::vector<int>::const_iterator itFind = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), iE);
	assert(itFind != vecProbToEdgeMapping.end());

	// Get the position (could directly use the iterator but illisible)
	size_t iVecProb = std::distance(vecProbToEdgeMapping.begin(), itFind);
	assert(iVecProb < vecScaling.size());
	return vecScaling[iVecProb];
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::setScalingFactorByVecPos(size_t iV, double aScalingFactor) {
	assert(iV < vecScaling.size());
	vecScaling[iV] = aScalingFactor;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
void EigenState<conditionalProbType>::setScalingFactorByEdgeId(size_t iE, double aScalingFactor) {
	// Look at the vec prob position
	std::vector<int>::iterator itFind = std::find(vecProbToEdgeMapping.begin(), vecProbToEdgeMapping.end(), iE);
	assert(itFind != vecProbToEdgeMapping.end());

	// Get the position (could directly use the iterator but illisible)
	size_t iVecProb = std::distance(vecProbToEdgeMapping.begin(), itFind);
	assert(iVecProb < vecScaling.size());
	vecScaling[iVecProb] = aScalingFactor;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
std::string EigenState<conditionalProbType>::toString() const {
	std::stringstream ss;

	ss << "ID =  " << id << std::endl;
	ss << "Number of probability vector allocated: ";
	ss << "[ Observed = " << (stateProb != NULL ? size() : 0)  << " ] " << std::endl;

	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
	if(stateProb != NULL) {
		//assert(observedStateProb->matrix.cols() == (int)vecProbToEdgeMapping.size());
		for(size_t iE=0; iE<N_VECTOR; ++iE) {
			ss << "U/S/UHat/SHat [" << iE << "] : " << stateProb->matrix.col(iE).format(CommaInitFmt) << std::endl;
		}

		for(size_t iE = 0; iE<vecProbToEdgeMapping.size(); ++iE) {
			ss << "Observed prob [" << iE << "] for edge " << vecProbToEdgeMapping[iE] << " - matId=" << stateProb->id << " : " << stateProb->matrix.col(iE+N_VECTOR).format(CommaInitFmt) << std::endl;
		}
	} else {
		assert(vecProbToEdgeMapping.empty());
	}


	return ss.str();
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
Likelihood::Monitor::ProbeState EigenState<conditionalProbType>::toProbe() const {
	Likelihood::Monitor::ProbeState probe;
	probe.u = getUnobservedStateProb();

	for(size_t i=0; i<size(); ++i) {
		probe.vecP.push_back(getObservedStateProb().col(i));
		probe.vecIdEdge.push_back(getVecProbToEdgeMapping()[i]);
		probe.scalingFactor.push_back(getScalingFactorByVecPos(i));
	}

	return probe;
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType> operator/( const EigenState<conditionalProbType> &state1 , const EigenState<conditionalProbType> &state2 ) {
	assert(state1.size() == state2.size());
	EigenState<conditionalProbType> newState(state1);

	//std::cout << " Element-wise division of id  = " << state1.id << " and id = " << state2.id << " into id = " << newState.id << std::endl;
	if(state1.stateProb != NULL) {
		newState.getStateProb() = state1.getStateProb().cwiseQuotient(state2.getStateProb());
	}

	return newState;
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType> abs( const EigenState<conditionalProbType> &state ) {
	EigenState<conditionalProbType> newState(state);

	//std::cout << " Abs of id  = " << state.id << " into  id = " << newState.id << std::endl;

	if(state.stateProb != NULL) {
		newState.getStateProb() = state.getStateProb().cwiseAbs();
	}

	return newState;
}


} /* namespace Optimized */
} /* namespace StateType */
} /* namespace Likelihood */

#endif // LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATE_HPP_

