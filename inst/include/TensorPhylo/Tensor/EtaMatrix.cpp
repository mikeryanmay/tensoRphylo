/*
 * EtaMatrix.cpp
 *
 *  Created on: Apr 2, 2020
 *      Author: meyerx
 */

#include "EtaMatrix.h"

#include <iostream>

namespace Tensor {

EtaMatrix::EtaMatrix() : etaStructure(ETA_DENSE) {
}

EtaMatrix::EtaMatrix(etaStructure_t aEtaStructure) : etaStructure(aEtaStructure) {
}

EtaMatrix::~EtaMatrix() {
}

void EtaMatrix::setStructureType(etaStructure_t aEtaStructure) {
	etaStructure = aEtaStructure;
}

etaStructure_t EtaMatrix::getStructureType() const {
	return etaStructure;
}

bool EtaMatrix::isQuasseMatrix(const Eigen::MatrixXd &aOriginalMatrix) const{
	assert(aOriginalMatrix.rows() == aOriginalMatrix.cols());

	double factor = 0;
	bool isQuasse = true;
	for(size_t i=0; i<(size_t)aOriginalMatrix.rows(); ++i) {
		if(i==0) {
			factor = aOriginalMatrix(i,i+1);
			isQuasse = isQuasse && fabs(aOriginalMatrix(i,i) + aOriginalMatrix(i,i+1)) < 1.e-10;
		} else if(i == (size_t)aOriginalMatrix.rows()-1) {
			isQuasse = isQuasse && fabs(aOriginalMatrix(i,i-1) + aOriginalMatrix(i,i)) < 1.e-10;
			isQuasse = isQuasse && fabs(aOriginalMatrix(i,i-1)-factor) < 1.e-10;
		} else {
			isQuasse = isQuasse && fabs(aOriginalMatrix(i,i-1) + aOriginalMatrix(i,i) + aOriginalMatrix(i,i+1)) < 1.e-10;
			isQuasse = isQuasse && fabs(aOriginalMatrix(i,i-1)-factor) < 1.e-10 && fabs(aOriginalMatrix(i,i)+2*factor) < 1.e-10;
		}
	}

	return isQuasse && aOriginalMatrix.rows() >= (1 << 4);
}


} /* namespace Tensor */
