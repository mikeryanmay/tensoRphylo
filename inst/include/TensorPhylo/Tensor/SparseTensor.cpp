/*
 * SparseTensor.cpp
 *
 *  Created on: March 9, 2020
 *      Author: meyerx
 */

#include "SparseTensor.h"

namespace Tensor {

SparseTensor::SparseTensor(size_t aN, const std::map< std::vector<unsigned>, double > &aEventMap) :
		BaseTensor(), N(aN), rbEventMap(aEventMap) {
	sparse = !rbEventMap.empty();
	init();
}

SparseTensor::~SparseTensor() {
}

tensor_t& SparseTensor::getTensor(double aT) {
	assert(false);
	static tensor_t dummyTensor;
	return dummyTensor;
}

void SparseTensor::doResetTensor() {
	dimensions.clear();
	sparseTensor.clear();
}

void SparseTensor::doDefineIfTimeDependent() {
	isTimeDependent = false;
}

void SparseTensor::doDefineDimensions() {

	dimensions.resize(3);
	if(rbEventMap.empty()) {
		dimensions[0] = 0;
		dimensions[1] = 0;
		dimensions[2] = 0;
	} else {
		dimensions[0] = N;
		dimensions[1] = N;
		dimensions[2] = N;
	}
}

void SparseTensor::doDefineTensor() {

	std::vector< std::vector < Eigen::Triplet<double> > > tripletsLists(N);

	// For each event, prepare for insertion
	for(cItRbEventMap_t itMap = rbEventMap.begin(); itMap != rbEventMap.end(); ++itMap) {
		assert(itMap->first.size() == dimensions.size());

		// Get x,y,z and p from an rbEventMap list
		size_t iParent = itMap->first[0];
		size_t iChildL =  itMap->first[1];
		size_t iChildR =  itMap->first[2];
		double prob = itMap->second;

		// Create an eigen triple
		Eigen::Triplet<double> aTriplet(iChildL, iChildR, prob);
		tripletsLists[iParent].push_back(aTriplet);
	}

	// Init all sparse matrix with the triple lists
	sparseTensor.resize(dimensions[0], Eigen::SparseMatrix< double, Eigen::RowMajor >(dimensions[1], dimensions[2]) );
	for(size_t iX=0; iX < dimensions[0]; ++iX) {
		sparseTensor[iX].setFromTriplets(tripletsLists[iX].begin(), tripletsLists[iX].end());
		sparseTensor[iX].makeCompressed();
	}
}

double SparseTensor::getNonZeroCoefficientRatio() const {

	size_t totalCoeff = dimensions[0]*dimensions[1]*dimensions[2];

	return (double)rbEventMap.size()/(double)totalCoeff;
}

bool SparseTensor::isSparse() const {
	return sparse;
}

sparseTensor_t& SparseTensor::getSparseTensor(double aT) {
	assert(sparse);
	return sparseTensor;
}

} /* namespace Tensor */
