/*
 * SimpleTensor.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: mrmay
 */

#include "SimpleTensor.h"

namespace Tensor {

SimpleTensor::SimpleTensor(std::vector<Eigen::MatrixXd> &aTensor) :
		BaseTensor(), originalTensor(aTensor) {
	sparse = false;
	init();
}

SimpleTensor::~SimpleTensor() {
}

tensor_t& SimpleTensor::getTensor(double aT) {
	return originalTensor;
}

void SimpleTensor::doResetTensor() {
	dimensions.clear();

	sparseTensor.clear();
}

void SimpleTensor::doDefineIfTimeDependent() {
	isTimeDependent = false;
}

void SimpleTensor::doDefineDimensions() {

	dimensions.resize(3);
	if(originalTensor.empty()) {
		dimensions[0] = 0;
		dimensions[1] = 0;
		dimensions[2] = 0;
	} else {
		dimensions[0] = originalTensor.size();
		dimensions[1] = originalTensor.at(0).rows();
		dimensions[2] = originalTensor.at(0).cols();
	}
}

void SimpleTensor::doDefineTensor() {

	sparse = !originalTensor.empty() && (getNonZeroCoefficientRatio() < NONZERO_COEFF_SPARSE_THRESHOLD);

	if(sparse) {
		for(size_t iX=0; iX < dimensions[0]; ++iX) {
			sparseTensor.push_back(originalTensor[iX].sparseView());
			sparseTensor[iX].makeCompressed();
		}
	}

}

double SimpleTensor::getNonZeroCoefficientRatio() const {

	size_t nonZeroCoeff = 0;
	size_t totalCoeff = dimensions[0]*dimensions[1]*dimensions[2];

	for(size_t iX=0; iX<dimensions[0]; ++iX) {
		for(size_t iY=0; iY<dimensions[1]; ++iY) {
			for(size_t iZ=0; iZ<dimensions[2]; ++iZ) {
				if(originalTensor[iX](iY, iZ) != 0.)  nonZeroCoeff++;
			}
		}
	}

	return (double)nonZeroCoeff/(double)totalCoeff;
}

bool SimpleTensor::isSparse() const {
	return sparse;
}

sparseTensor_t& SimpleTensor::getSparseTensor(double aT) {
	assert(sparse);
	return sparseTensor;
}

} /* namespace Tensor */
