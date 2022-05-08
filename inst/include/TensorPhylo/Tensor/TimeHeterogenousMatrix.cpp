/*
 * TimeHeterogenousMatrix.cpp
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#include "TimeHeterogenousMatrix.h"

#include <iostream>

namespace Tensor {

TimeHeterogenousMatrix::TimeHeterogenousMatrix(
		const std::vector<double> &aTimes,
		const std::vector< Eigen::MatrixXd > &aParameters) :
		BaseTensor(), EtaMatrix(), times(aTimes), parameters(aParameters) {
	sparse = false;

	assert(times.size()+1 == parameters.size());

	init();
}

TimeHeterogenousMatrix::~TimeHeterogenousMatrix() {
}

tensor_t& TimeHeterogenousMatrix::getTensor(double aT) {
	size_t epochIndex = getEpochIndex(aT, times);
	assert(epochIndex < vecTensors.size());

	return vecTensors[epochIndex];
}

sparseTensor_t& TimeHeterogenousMatrix::getSparseTensor(double aT) {
	assert(sparse);

	size_t epochIndex = getEpochIndex(aT, times);
	assert(epochIndex < vecSparseTensor.size());

	return vecSparseTensor[epochIndex];
}

void TimeHeterogenousMatrix::doResetTensor() {
	dimensions.clear();
	vecTensors.clear();

	vecSparseTensor.clear();
}

void TimeHeterogenousMatrix::doDefineIfTimeDependent() {
	isTimeDependent = true;
}

void TimeHeterogenousMatrix::doDefineDimensions() {
	dimensions.push_back(parameters.front().rows());
	dimensions.push_back(parameters.front().cols());
}

void TimeHeterogenousMatrix::doDefineTensor() {

	for(size_t iP=0; iP<parameters.size(); ++iP) {
		tensor_t tensor;
		tensor.push_back(parameters[iP]);
		tensor.back().diagonal().setZero();
		Eigen::VectorXd tmpVector = tensor.back().rowwise().sum();
		tensor.back().diagonal() = -1.0 * tmpVector;
		vecTensors.push_back(tensor);
	}

	sparse = getNonZeroCoefficientRatio() < NONZERO_COEFF_SPARSE_THRESHOLD;

	if(sparse) {
		bool isQuasse = true;
		for(size_t iP=0; iP<parameters.size(); ++iP) {
			{ // Create sparse view of matrix iP
			sparseTensor_t tmpSpTensor(1, vecTensors[iP][0].sparseView());
			tmpSpTensor.front().makeCompressed();
			vecSparseTensor.push_back(tmpSpTensor);
			}
			isQuasse = isQuasse && isQuasseMatrix(vecTensors[iP].front());
		}


		if(isQuasse) {
			setStructureType(ETA_QUASSE);
		} else {
			setStructureType(ETA_SPARSE);
			//setStructureType(ETA_DENSE);
		}

	} else {
		setStructureType(ETA_DENSE);
	}
}

double TimeHeterogenousMatrix::getNonZeroCoefficientRatio() const {

	size_t totalCoeff = dimensions[0]*dimensions[1];

	bool isSparse = true;

	double avgNonZeroCoefficientRatio = 0.;
	for(size_t iV=0; iV < vecTensors.size(); ++iV) {
		size_t nonZeroCoeff = 0;
		for(size_t iX=0; iX < dimensions[0]; ++iX) {
			for(size_t iY=0; iY < dimensions[1]; ++iY) {
				if(vecTensors[iV].front()(iX, iY) != 0.) nonZeroCoeff++;
			}
		}
		avgNonZeroCoefficientRatio += (double)nonZeroCoeff/(double)totalCoeff;
		isSparse = isSparse && (((double)nonZeroCoeff/(double)totalCoeff) < NONZERO_COEFF_SPARSE_THRESHOLD);
	}

	avgNonZeroCoefficientRatio /= vecTensors.size();
	if(isSparse && avgNonZeroCoefficientRatio < NONZERO_COEFF_SPARSE_THRESHOLD) {
		return avgNonZeroCoefficientRatio;
	} else if(isSparse && avgNonZeroCoefficientRatio >= NONZERO_COEFF_SPARSE_THRESHOLD) {
		return NONZERO_COEFF_SPARSE_THRESHOLD/2.0; // hack
	} else if(!isSparse && avgNonZeroCoefficientRatio < NONZERO_COEFF_SPARSE_THRESHOLD) {
		return 1.0; // Hack
	} else {
		return avgNonZeroCoefficientRatio;
	}

}


bool TimeHeterogenousMatrix::isSparse() const {
	return sparse;
}



} /* namespace Tensor */
