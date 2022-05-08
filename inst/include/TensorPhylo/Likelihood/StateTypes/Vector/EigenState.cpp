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
namespace Vector {

size_t EigenState::idSeq = 0;

EigenState::EigenState() : id(idSeq++), edgeMapping(-1), scaling(0.), probVec(NULL) {
}

EigenState::EigenState(const EigenState &aEigenState) :
		id(idSeq++), edgeMapping(aEigenState.edgeMapping), scaling(aEigenState.scaling), probVec(NULL) {

	if(aEigenState.probVec != NULL) {
		initMemory();
		probVec->vector = aEigenState.probVec->vector;
	} else {
		probVec = NULL;
	}

}

EigenState::~EigenState() {
	//std::cout << " Delete ID = " << id << std::endl;

	if(probVec != NULL) {
		Utils::MemoryPool::eigenCPU().threadSafeFreeVector(probVec);
		probVec = NULL;
	}

}

void EigenState::setEdgeMapping(long int aEdgeId) {
	edgeMapping = aEdgeId;
}

long int EigenState::getEdgeMapping() const {
	return edgeMapping;
}

void EigenState::resize() {

	// Do nothing instead?
	if(probVec == NULL) {
		allocateVecProb();
	}
}

void EigenState::allocateVecProb() {
	assert(probVec == NULL);
	if(probVec == NULL) {
		probVec = Utils::MemoryPool::eigenCPU().threadSafeAllocateVector();
	}
}

void EigenState::removeVecProb() {
	assert(probVec != NULL);
	if(probVec != NULL) {
		Utils::MemoryPool::eigenCPU().threadSafeFreeVector(probVec);
		probVec = NULL;
	}
}


Eigen::VectorXd& EigenState::getStateProb() {
	if(probVec == NULL) {
		probVec = Utils::MemoryPool::eigenCPU().threadSafeAllocateVector();
	}
	return probVec->vector;
}

const Eigen::VectorXd& EigenState::getStateProb() const {
	assert(probVec != NULL);
	return probVec->vector;
}

void EigenState::setScaling(double aScaling) {
	scaling = aScaling;
}

double EigenState::getScaling() const {
	return scaling;
}

size_t EigenState::size() const {
	return 1;
}

double EigenState::defineNormInf() const {
	//std::cout << " Norm inf of id  = " << id << std::endl;

    double absMax = (double)probVec->vector.lpNorm<Eigen::Infinity>();
    return absMax;
}

void EigenState::addMult(double factor, const EigenState &otherState ) {
	assert(otherState.size() == size());
	probVec->vector += factor*otherState.probVec->vector;
}

void EigenState::initMult(double factor, const EigenState &otherState ) {
	resize();

	probVec->vector = factor*otherState.probVec->vector;
	scaling = otherState.scaling;
	edgeMapping = otherState.edgeMapping;
}

void EigenState::odeIntRelativeError(double m_eps_abs, double m_eps_rel, double m_a_x, double m_a_dxdt, const EigenState &aState1, const EigenState &aState2) {
	assert(aState1.size() == size() && aState2.size() == size());

	Eigen::VectorXd denom = (m_eps_abs + (m_eps_rel * ( m_a_x * aState1.probVec->vector.cwiseAbs() + m_a_dxdt * aState2.probVec->vector.cwiseAbs())).array()).matrix();
	probVec->vector = probVec->vector.cwiseAbs().cwiseQuotient(denom);
}

EigenState& EigenState::operator+=( const double &val ) {
	//std::cout << " ADD val = " << val << " to id = " << id << std::endl;
	probVec->vector = (probVec->vector.array() + val).matrix();

	return *this;
}

EigenState& EigenState::operator+=( const EigenState &otherState ) {

	//std::cout << " ADD id = " << probVec.id << " to id = " << id << std::endl;
	assert(otherState.size() == size());
	probVec->vector += otherState.probVec->vector;
	return *this;
}

EigenState& EigenState::operator*=( const double a ) {
	//std::cout << " Multiply id = " << id << " by = " << a << std::endl;

	probVec->vector *= a;

	return *this;
}


EigenState& EigenState::operator=(const EigenState &otherState) {
	//std::cout << " Affect = " << probVec.id << " to id = " << id << std::endl;
	if(otherState.probVec != NULL) {
		resize();

		probVec->vector = otherState.probVec->vector;
		edgeMapping = otherState.edgeMapping;
		scaling = otherState.scaling;
	}

	return *this;
}

void EigenState::initMemory() {
	assert(probVec == NULL);
	probVec = Utils::MemoryPool::eigenCPU().threadSafeAllocateVector();
}


void EigenState::roundNegativeProbabilityToZero() {

	if((probVec->vector.array() < 0.).any()) {
		for(size_t iS=0; iS<(size_t)probVec->vector.size(); ++iS) {
			if(probVec->vector(iS) < 0.) {
				probVec->vector(iS) = 0.;
			}
		}
	}

}

void EigenState::rescaleProbabilities() {
	scaling = rescaleProbabilityVector(probVec->vector);
}


void EigenState::copyScalingFactors(const EigenState &otherState) {
	scaling = otherState.scaling;
}

std::string EigenState::toString() const {
	std::stringstream ss;

	ss << "ID =  " << id << std::endl;

	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
	if(probVec) {
		ss << "Prob - vecId=" << probVec->id << " : " << probVec->vector.format(CommaInitFmt) << std::endl;
	}

	return ss.str();
}

Likelihood::Monitor::ProbeState EigenState::toProbe() const {
	Likelihood::Monitor::ProbeState probe;

	probe.vecP.push_back(probVec->vector);
	probe.vecIdEdge.push_back(edgeMapping);
	probe.scalingFactor.push_back(scaling);

	return probe;
}


EigenState operator/( const EigenState &otherState1 , const EigenState &otherState2 ) {
	assert(otherState1.size() == otherState2.size());
	EigenState newState(otherState1);

	//std::cout << " Element-wise division of id  = " << probVec1.id << " and id = " << probVec2.id << " into id = " << newState.id << std::endl;

	newState.probVec->vector = otherState1.probVec->vector.cwiseQuotient(otherState2.probVec->vector);

	return newState;
}

EigenState abs( const EigenState &otherState ) {
	EigenState newState(otherState);

	//std::cout << " Abs of id  = " << probVec.id << " into  id = " << newState.id << std::endl;

	newState.probVec->vector = otherState.probVec->vector.cwiseAbs();

	return newState;
}

} /* namespace Default */
} /* namespace StateType */
} /* namespace Likelihood */

