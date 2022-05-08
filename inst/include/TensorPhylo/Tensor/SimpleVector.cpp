/*
 * SimpleVector.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#include "SimpleVector.h"

namespace Tensor {

SimpleVector::SimpleVector(Eigen::VectorXd &aVector) :
		BaseTensor(), originalVector(aVector) {
	sparse = false;
	init();
}

SimpleVector::~SimpleVector() {
}

tensor_t& SimpleVector::getTensor(double aT) {
	return tensor;
}


void SimpleVector::doResetTensor() {
	dimensions.clear();
	tensor.clear();
}

void SimpleVector::doDefineIfTimeDependent() {
	isTimeDependent = false;
}

void SimpleVector::doDefineDimensions() {
	dimensions.assign(1, originalVector.size());
}

void SimpleVector::doDefineTensor() {
	tensor.push_back(originalVector);

	sparse = false; //getNonZeroCoefficientRatio() < NONZERO_COEFF_SPARSE_THRESHOLD;
}

double SimpleVector::getNonZeroCoefficientRatio() const {

	size_t nonZeroCoeff = 0;
	size_t totalCoeff = dimensions[0];

	for(size_t iX=0; iX<dimensions[0]; ++iX) {
		if(originalVector(iX) != 0.)  nonZeroCoeff++;
	}

	return (double)nonZeroCoeff/(double)totalCoeff;
}


bool SimpleVector::isSparse() const {
	return sparse;
}



} /* namespace Tensor */
