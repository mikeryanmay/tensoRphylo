/*
 * TimeHeterogenousVector.cpp
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#include "TimeHeterogenousVector.h"

#include <iostream>

namespace Tensor {

TimeHeterogenousVector::TimeHeterogenousVector(
		const std::vector<double> &aTimes,
		const std::vector< Eigen::VectorXd > &aParameters) :
		BaseTensor(), times(aTimes), parameters(aParameters) {
	sparse = false;

	assert(times.size()+1 == parameters.size());

	init();
}

TimeHeterogenousVector::~TimeHeterogenousVector() {
}

tensor_t& TimeHeterogenousVector::getTensor(double aT) {
	size_t epochIndex = getEpochIndex(aT, times);
	assert(epochIndex < vecTensors.size());

	/*if(times.size() > 0) {
		std::cout << "Time : " << std::scientific << aT << " -> " << times[0] << " -> " << aT-times[0] << " :: "<< std::fixed << epochIndex << std::endl;
	}*/

	return vecTensors[epochIndex];
}

void TimeHeterogenousVector::doResetTensor() {
	dimensions.clear();
	vecTensors.clear();
}

void TimeHeterogenousVector::doDefineIfTimeDependent() {
	isTimeDependent = true;
}

void TimeHeterogenousVector::doDefineDimensions() {
	dimensions.assign(1, parameters.front().size());
}

void TimeHeterogenousVector::doDefineTensor() {

	vecTensors.clear();

	for(size_t iP=0; iP<parameters.size(); ++iP) {
		tensor_t tensor;
		tensor.push_back(parameters[iP]);
		vecTensors.push_back(tensor);
	}

	sparse = false; //getNonZeroCoefficientRatio() < NONZERO_COEFF_SPARSE_THRESHOLD;
}

double TimeHeterogenousVector::getNonZeroCoefficientRatio() const {

	size_t nonZeroCoeff = 0;
	size_t totalCoeff = dimensions[0]*vecTensors.size();

	bool isSparse = true;

	for(size_t iV=0; iV<vecTensors.size(); ++iV) {
		for(size_t iX=0; iX<dimensions[0]; ++iX) {
			if(vecTensors[iV].front()(iX) != 0.) nonZeroCoeff++;
		}
		isSparse = isSparse && ((double)nonZeroCoeff/(double)totalCoeff);
	}

	return isSparse;
}


bool TimeHeterogenousVector::isSparse() const {
	return sparse;
}



} /* namespace Tensor */
