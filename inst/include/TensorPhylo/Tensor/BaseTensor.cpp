/*
 * BaseTensor.cpp
 *
 *  Created on: Aug 28, 2019
 *      Author: xaviermeyer
 */

#include "BaseTensor.h"

namespace Tensor {

const double BaseTensor::NONZERO_COEFF_SPARSE_THRESHOLD = 0.5;

BaseTensor::BaseTensor() : isTimeDependent(false) {
}

BaseTensor::~BaseTensor() {
}

bool BaseTensor::isContantThroughTime() const {
	return !isTimeDependent;
}

std::vector<size_t> BaseTensor::getDimensions() {
	return dimensions;
}

sparseTensor_t& BaseTensor::getSparseTensor(double aT) {
	assert(false && "Not implemented for this specific type of tensors.");
	static sparseTensor_t dummy;
	return dummy;
}

void BaseTensor::update() {
	doResetTensor();
	init();
}


void BaseTensor::init() {
	doDefineIfTimeDependent();
	doDefineDimensions();
	doDefineTensor();
}

size_t BaseTensor::getEpochIndex(double currentTime, const std::vector<double> &aTimes) {

	//assert(!aTimes.empty() && "SimplexXXXXX tensor class should be used for time-homogeneous vectors.");

	for(size_t iT=0; iT<aTimes.size(); ++iT) {
		if(currentTime < aTimes[iT]) {
			return iT;
		}
	}
	return aTimes.size();
}

} /* namespace Tensor */
