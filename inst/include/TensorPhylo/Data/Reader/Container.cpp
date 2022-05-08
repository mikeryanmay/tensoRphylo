/*
 * Container.cpp
 *
 *  Created on: Mar 10, 2020
 *      Author: meyerx
 */

#include "Container.h"

#include <cassert>
#include <iostream>

#include "Utils/Parallel/Manager.h"

namespace Phylogeny {
namespace Data {

Container::Container() {
}

Container::~Container() {
}

bool Container::isReady() const {
	return !mapTaxaToId.empty() &&
			!mapTaxaToProbs.empty() &&
			!mapIdToTaxa.empty() &&
			//mapTaxaToId.size() == mapTaxaToProbs.size() &&
			mapIdToTaxa.size() == mapTaxaToId.size();
}

void Container::setTaxaToIdMap(const std::vector<std::string> &aTaxaVec) {

	mapIdToTaxa.clear();
	mapTaxaToId.clear();

	for(size_t iT=0; iT<aTaxaVec.size(); ++iT) {
		mapIdToTaxa[iT] = aTaxaVec[iT];
		mapTaxaToId[aTaxaVec[iT]] = iT;
	}
}


void Container::setTaxaMaps(const std::map<std::string, std::size_t> &aTaxaToIdMap, const std::map<std::size_t, std::string> &aIdToTaxaMap) {
	mapIdToTaxa = aIdToTaxaMap;
	mapTaxaToId = aTaxaToIdMap;
}

void Container::registerTaxaData(const std::string &taxaName, const std::vector<double> &probs) {

	Eigen::VectorXd eigenProbs(probs.size());

	for(size_t iP=0; iP<probs.size(); ++iP) {
		eigenProbs(iP) = probs[iP];
	}

	mapTaxaToProbs[taxaName] = eigenProbs;
}

void Container::registerTaxaData(size_t nStates, const std::string &taxaName, const std::vector<int> &states) {

	Eigen::VectorXd eigenProbs(nStates);

	assert(!states.empty());
	if(states.size() == 1 && states.front() == -1) {
		eigenProbs.setOnes();
	} else {
		eigenProbs.setZero();
		for(size_t iS=0; iS<states.size(); ++iS) {
			assert(states[iS] >= 0 && states[iS] < (long int)nStates);
			eigenProbs(states[iS]) = 1.0;
		}
	}
	mapTaxaToProbs[taxaName] = eigenProbs;
}

size_t Container::getPositionForTaxaId(int idTaxa) const {
	assert(isReady());

	// this is ugly but more efficient to get the states
	for(size_t iF=0; iF<vecProbCache.size(); iF++) {
		if(vecProbCache[iF].idTaxa == idTaxa) {
			return iF;
		}
	}
	vecProbCache.reserve(mapIdToTaxa.size());

	// we didn't find the data
	stateCache_t stateCache;
	stateCache.idTaxa = idTaxa;

	std::map< size_t, std::string >::const_iterator itFindTaxa = mapIdToTaxa.find(idTaxa);
	assert(itFindTaxa != mapIdToTaxa.end() && "A taxa was not found was not found.");

	std::string taxaLabel = itFindTaxa->second;
	std::map< std::string, Eigen::VectorXd  >::const_iterator itFindProbs = mapTaxaToProbs.find(taxaLabel);
	assert(itFindProbs != mapTaxaToProbs.end() && "State probabilities for a taxa was not found.");

	stateCache.statesProb = itFindProbs->second;

	size_t iSC = vecProbCache.size();
	vecProbCache.push_back(stateCache);

	return iSC;
}

const Eigen::VectorXd& Container::getProbForTaxaId(int idTaxa) const {
	assert(isReady());
	size_t iPos = getPositionForTaxaId(idTaxa);
	return vecProbCache[iPos].statesProb ;
}

const Eigen::VectorXd& Container::getProbForTaxaIdThreadSafe(int idTaxa) const {
#ifdef _OPENMP
	if(!omp_in_parallel()) {
		return getProbForTaxaId(idTaxa);
	} else {
		size_t iPos = 0;
		#pragma omp critical (data_probability_cache)
		{
			iPos = getPositionForTaxaId(idTaxa);
		}
		return vecProbCache[iPos].statesProb ;
	}
#else
	return getProbForTaxaId(idTaxa);
#endif
}

const Eigen::VectorXd&  Container::getProbForTaxaLabel(const std::string &taxaLabel) const {
	assert(isReady());

	std::map< std::string, Eigen::VectorXd  >::const_iterator itFind = mapTaxaToProbs.find(taxaLabel);
	assert(itFind != mapTaxaToProbs.end() && "State probabilities for a taxa was not found.");

	return itFind->second;
}

} /* namespace Data */
} /* namespace Phylogeny */
