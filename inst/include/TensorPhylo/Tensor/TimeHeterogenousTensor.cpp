/*
 * TimeHeterogenousTensor.cpp
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#include "TimeHeterogenousTensor.h"

#include <iostream>

namespace Tensor {

TimeHeterogenousTensor::TimeHeterogenousTensor(
		size_t aN,
		const std::vector<double> &aTimes,
		const std::vector< rbEventMap_t > &aParameters) :
		BaseTensor(), N(aN), times(aTimes), parameters(aParameters) {

	empty = (times.empty() && parameters.empty());
	assert(((times.size()+1) == parameters.size()) || empty);
	sparse = !empty;

	init();
}

TimeHeterogenousTensor::~TimeHeterogenousTensor() {
}

tensor_t& TimeHeterogenousTensor::getTensor(double aT) {
	assert(false);
	static tensor_t dummyTensor;
	return dummyTensor;
}

void TimeHeterogenousTensor::doDefineIfTimeDependent() {
	isTimeDependent = true;
}


sparseTensor_t& TimeHeterogenousTensor::getSparseTensor(double aT) {
	assert(!empty);
	size_t epochIndex = getEpochIndex(aT, times);
	assert(epochIndex < vecSparseTensors.size());

	return vecSparseTensors[epochIndex];
}

void TimeHeterogenousTensor::doResetTensor() {
	dimensions.clear();
	vecSparseTensors.clear();
}

void TimeHeterogenousTensor::doDefineDimensions() {

	dimensions.resize(3);
	if(parameters.empty()) {
		dimensions[0] = 0;
		dimensions[1] = 0;
		dimensions[2] = 0;
	} else {
		dimensions[0] = N;
		dimensions[1] = N;
		dimensions[2] = N;
	}
}

void TimeHeterogenousTensor::doDefineTensor() {

	vecSparseTensors.resize(parameters.size());

	for(size_t iP=0; iP < parameters.size(); ++iP) {

		std::vector< std::vector < Eigen::Triplet<double> > > tripletsLists(N);

		// For each event, prepare for insertion
		for(cItRbEventMap_t itMap = parameters[iP].begin(); itMap != parameters[iP].end(); ++itMap) {
			assert(itMap->first.size() == dimensions.size());

			// Get x,y,z and p from an rbEventMap list
			size_t iParent = itMap->first[0];
			size_t iChildL = itMap->first[1];
			size_t iChildR = itMap->first[2];
			double prob = itMap->second;

			// Create an eigen triple
			Eigen::Triplet<double> aTriplet(iChildL, iChildR, prob);
			tripletsLists[iParent].push_back(aTriplet);
		}

		// Init all sparse matrix with the triple lists
		vecSparseTensors[iP].resize(dimensions[0], Eigen::SparseMatrix< double, Eigen::RowMajor >(dimensions[1], dimensions[2]) );
		for(size_t iX=0; iX < dimensions[0]; ++iX) {
			vecSparseTensors[iP][iX].setFromTriplets(tripletsLists[iX].begin(), tripletsLists[iX].end());
			vecSparseTensors[iP][iX].makeCompressed();
		}

	}
}

double TimeHeterogenousTensor::getNonZeroCoefficientRatio() const {

	size_t totalCoeff = dimensions[0]*dimensions[1]*dimensions[2];

	size_t sum = 0;
	for(size_t iP=0; iP<parameters.size(); ++iP) {
		sum += parameters[iP].size();
	}
	double mean = (double)sum/(double)parameters.size();

	return mean/(double)totalCoeff;
}

bool TimeHeterogenousTensor::isSparse() const {
	return sparse;
}


} /* namespace Tensor */
