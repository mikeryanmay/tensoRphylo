/*
 * SimpleMatrix.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: mrmay
 */

#include "SimpleMatrix.h"

#include <iostream>
namespace Tensor {

SimpleMatrix::SimpleMatrix(Eigen::MatrixXd &aMatrix) :
		BaseTensor(), EtaMatrix(), originalMatrix(aMatrix) {
	sparse = false;
	init();
}

SimpleMatrix::~SimpleMatrix() {
}

tensor_t& SimpleMatrix::getTensor(double aT) {
	return tensor;
}

void SimpleMatrix::doResetTensor() {
	dimensions.clear();
	tensor.clear();

	sparseTensor.clear();
}

void SimpleMatrix::doDefineIfTimeDependent() {
	isTimeDependent = false;
}

void SimpleMatrix::doDefineDimensions() {
	dimensions.push_back(originalMatrix.rows());
	dimensions.push_back(originalMatrix.cols());
}

void SimpleMatrix::doDefineTensor() {
	tensor.push_back(originalMatrix);
	originalMatrix.diagonal().setZero();
	Eigen::VectorXd tmpVector = originalMatrix.rowwise().sum();
	originalMatrix.diagonal() = -1.0 * tmpVector;


	sparse = getNonZeroCoefficientRatio() < NONZERO_COEFF_SPARSE_THRESHOLD;

	bool isQuasse = sparse && isQuasseMatrix(originalMatrix);

	if(sparse) {
		sparseTensor.push_back(tensor[0].sparseView());
		sparseTensor[0].makeCompressed();
	}

	if(sparse && !isQuasse) {
		setStructureType(ETA_SPARSE);
	} else if(sparse && isQuasse) {
		setStructureType(ETA_QUASSE);
	} else {
		setStructureType(ETA_DENSE);
	}

}

double SimpleMatrix::getNonZeroCoefficientRatio() const {

	size_t nonZeroCoeff = 0;
	size_t totalCoeff = dimensions[0]*dimensions[1];

	for(size_t iX=0; iX<dimensions[0]; ++iX) {
		for(size_t iY=0; iY<dimensions[1]; ++iY) {
			if(originalMatrix(iX, iY) != 0.)  nonZeroCoeff++;
		}
	}

	return (double)nonZeroCoeff/(double)totalCoeff;
}

bool SimpleMatrix::isSparse() const {
	return sparse;
}

sparseTensor_t& SimpleMatrix::getSparseTensor(double aT) {
	assert(sparse);
	return sparseTensor;
}

} /* namespace Tensor */
