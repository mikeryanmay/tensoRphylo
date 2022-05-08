/*
 * EigenState.cpp
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#include "EigenState.h"

#include "Utils/MemoryPool/EigenCPU.h"
#include "../Utils.h"

using Likelihood::StateType::rescaleProbabilityVector;

namespace Likelihood {
namespace StateType {
namespace Matrix {

size_t EigenState::idSeq = 0;

EigenState::EigenState() : id(idSeq++){
}

EigenState::EigenState(const EigenState &aEigenState) :
		id(idSeq++) {

	initMemory(aEigenState);
}

EigenState::~EigenState() {
	//std::cout << " Delete ID = " << id << std::endl;

	for(size_t i=0; i<probVec.size(); ++i) {
		Utils::MemoryPool::eigenCPU().threadSafeFreeVector(probVec[i]);
	}
	probVec.clear();

}

void EigenState::resize(size_t aSize) {

	size_t initSize = size();

	if(aSize < initSize) {
		for(size_t i=aSize; i<size(); ++i) {
			Utils::MemoryPool::eigenCPU().threadSafeFreeVector(probVec[i]);
		}
	}

	probVec.resize(aSize);

	for(size_t i=initSize; i<aSize; ++i) {
		probVec[i] = Utils::MemoryPool::eigenCPU().threadSafeAllocateVector();
	}
}

Eigen::VectorXd& EigenState::getStateProb(size_t iPos) {
	if(iPos >= size()) {
		resize(iPos+1);
	}
	return probVec[iPos]->vector;
}

const Eigen::VectorXd& EigenState::getStateProb(size_t iPos) const {
	assert(iPos < size());
	return probVec[iPos]->vector;
}

size_t EigenState::size() const {
	return probVec.size();
}

double EigenState::defineNormInf() const {
	//std::cout << " Norm inf of id  = " << id << std::endl;

    double absMax = 0.;
    for(size_t iV=0; iV<size(); ++iV) {
    	absMax = std::max(absMax, (double)probVec[iV]->vector.lpNorm<Eigen::Infinity>());
    }
    return absMax;
}

void EigenState::addMult(double factor, const EigenState &otherState ) {
	assert(otherState.size() == size());
    for(size_t iV=0; iV<size(); ++iV) {
    	probVec[iV]->vector += factor*otherState.probVec[iV]->vector;
    }
}

void EigenState::initMult(double factor, const EigenState &otherState ) {
	resize(otherState.size());

	for(size_t iV=0; iV<size(); ++iV) {
		probVec[iV]->vector = factor*otherState.probVec[iV]->vector;
	}

}

void EigenState::odeIntRelativeError(double m_eps_abs, double m_eps_rel, double m_a_x, double m_a_dxdt, const EigenState &aState1, const EigenState &aState2) {
	assert(aState1.size() == size() && aState2.size() == size());

	for(size_t iV=0; iV<size(); ++iV) {
		Eigen::VectorXd denom = (m_eps_abs + (m_eps_rel * ( m_a_x * aState1.probVec[iV]->vector.cwiseAbs() + m_a_dxdt * aState2.probVec[iV]->vector.cwiseAbs())).array()).matrix();
		probVec[iV]->vector = probVec[iV]->vector.cwiseAbs().cwiseQuotient(denom);
	}
}

EigenState& EigenState::operator+=( const double &val ) {
	//std::cout << " ADD val = " << val << " to id = " << id << std::endl;
	for(size_t iV=0; iV<size(); ++iV) {
		probVec[iV]->vector = (probVec[iV]->vector.array() + val).matrix();
	}

	return *this;
}

EigenState& EigenState::operator+=( const EigenState &otherState ) {

	//std::cout << " ADD id = " << probVec.id << " to id = " << id << std::endl;
	assert(otherState.size() == size());
	for(size_t iV=0; iV<size(); ++iV) {
		probVec[iV]->vector += otherState.probVec[iV]->vector;
	}
	return *this;
}

EigenState& EigenState::operator*=( const double a ) {
	//std::cout << " Multiply id = " << id << " by = " << a << std::endl;

	for(size_t iV=0; iV<size(); ++iV) {
		probVec[iV]->vector *= a;
	}

	return *this;
}


EigenState& EigenState::operator=(const EigenState &otherState) {
	//std::cout << " Affect = " << probVec.id << " to id = " << id << std::endl;
	if(!otherState.probVec.empty()) {
		resize(otherState.size());

		for(size_t iV=0; iV<size(); ++iV) {
			probVec[iV]->vector = otherState.probVec[iV]->vector;
		}
	}

	return *this;
}

void EigenState::initMemory(const EigenState &aEigenState) {

	probVec.resize(aEigenState.probVec.size());
	for(size_t i=0; i<probVec.size(); ++i) {
		probVec[i] = Utils::MemoryPool::eigenCPU().threadSafeAllocateVector();
		probVec[i]->vector = aEigenState.probVec[i]->vector;
	}
}


void EigenState::roundNegativeProbabilityToZero() {

	for(size_t iV=0; iV<size(); ++iV) {
		if((probVec[iV]->vector.array() < 0.).any()) {
			for(size_t iS=0; iS<(size_t)probVec[iV]->vector.size(); ++iS) {
				if(probVec[iV]->vector(iS) < 0.) {
					probVec[iV]->vector(iS) = 0.;
				}
			}
		}
	}

}

std::string EigenState::toString() const {
	std::stringstream ss;

	ss << "ID =  " << id << std::endl;

	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
	for(size_t iV=0; iV<size(); ++iV) {
		ss << "Prob - vecId=" << probVec[iV]->id << " : " << probVec[iV]->vector.format(CommaInitFmt) << std::endl;
	}

	return ss.str();
}

Likelihood::Monitor::ProbeState EigenState::toProbe() const {
	Likelihood::Monitor::ProbeState probe;

	for(size_t iV=0; iV<size(); ++iV) {
		probe.vecP.push_back(probVec[iV]->vector);
		probe.vecIdEdge.push_back(-1);
		probe.scalingFactor.push_back(0.);
	}

	return probe;
}


EigenState operator/( const EigenState &otherState1 , const EigenState &otherState2 ) {
	assert(otherState1.size() == otherState2.size());
	EigenState newState(otherState1);

	//std::cout << " Element-wise division of id  = " << probVec1.id << " and id = " << probVec2.id << " into id = " << newState.id << std::endl;

	for(size_t iV=0; iV<otherState1.size(); ++iV) {
		newState.probVec[iV]->vector = otherState1.probVec[iV]->vector.cwiseQuotient(otherState2.probVec[iV]->vector);
	}

	return newState;
}

EigenState abs( const EigenState &otherState ) {
	EigenState newState(otherState);

	//std::cout << " Abs of id  = " << probVec.id << " into  id = " << newState.id << std::endl;
	for(size_t iV=0; iV<otherState.size(); ++iV) {
		newState.probVec[iV]->vector = otherState.probVec[iV]->vector.cwiseAbs();
	}

	return newState;
}

} /* namespace Matrix */
} /* namespace StateType */
} /* namespace Likelihood */

