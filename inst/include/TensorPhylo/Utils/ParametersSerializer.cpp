/*
 * ParametersSerializer.cpp
 *
 *  Created on: Mar 15, 2020
 *      Author: meyerx
 */

#include "ParametersSerializer.h"

#include <iostream>
#include <sstream>
#include <fstream>

namespace Utils {
namespace Serializer {
namespace Parameters {

std::string toString(const std::vector<double> &aVector) {
	std::stringstream ss;
	ss << aVector.size() << std::endl;
	for(size_t i=0; i<aVector.size(); ++i) {
		ss << aVector[i];
		if(i != aVector.size()-1) ss << " ";
	}
	if(!aVector.empty()) ss << std::endl;

	return ss.str();
}

std::string toString(const std::vector<Eigen::VectorXd> &aVector2) {
	std::stringstream ss;
	ss << aVector2.size() << std::endl;
	for(size_t i=0; i<aVector2.size(); ++i) {
		ss << toString(aVector2[i]);
	}
	return ss.str();
}

std::string toString(const std::vector<Eigen::MatrixXd> &aTensor) {
	std::stringstream ss;
	ss << aTensor.size();
	if(aTensor.size() == 0) {
		ss << std::endl;
		return ss.str();
	}

	ss << " " <<  aTensor[0].rows() <<  " " << aTensor[0].cols() << std::endl;
	for(size_t k = 0; k<aTensor.size(); ++k) {
		for(size_t i = 0; i < (size_t)aTensor[k].rows(); ++i) {
			for(size_t j = 0; j < (size_t)aTensor[k].cols(); ++j) {
				ss << aTensor[k](i,j);
				if(j != (size_t)aTensor[k].cols()-1) ss << " ";
			}
			ss << std::endl;
		}
	}
	return ss.str();
}

std::string toString(const Eigen::VectorXd &aVector) {
	std::stringstream ss;
	ss << aVector.size() << std::endl;
	for(size_t i=0; i<(size_t)aVector.size(); ++i) {
		ss << aVector(i);
		if(i != (size_t)aVector.size()-1) ss << " ";
	}
	ss << std::endl;

	return ss.str();
}

std::string toString(const Eigen::MatrixXd &aMatrix) {
	std::stringstream ss;
	ss << aMatrix.rows() << " "  << aMatrix.cols() << std::endl;
	for(size_t i=0; i<(size_t)aMatrix.rows(); ++i) {
		for(size_t j=0; j<(size_t)aMatrix.cols(); ++j) {
			ss << aMatrix(i,j);
			if(j != (size_t)aMatrix.cols()-1) ss << " ";
		}
		ss << std::endl;
	}
	return ss.str();
}


std::string toString(const std::vector< Tensor::rbEventMap_t > &aVecSparseTensor) {
	std::stringstream ss;
	ss << aVecSparseTensor.size() << std::endl;
	for(size_t i=0; i<aVecSparseTensor.size(); ++i) {
		ss << aVecSparseTensor[i].size() << std::endl;
		for(Tensor::rbEventMap_t::const_iterator it=aVecSparseTensor[i].begin(); it != aVecSparseTensor[i].end(); ++it) {
			assert(it->first.size() == 3);
			ss << it->first[0] << " " << it->first[1] << " " << it->first[2] << " " << it->second << std::endl;
		}
	}
	return ss.str();
}


void fromFile(std::ifstream &myFile,  std::vector<double> &aVector) {

	size_t size = 0;
	myFile >> size;
	aVector.clear();

	for(size_t i=0; i<size; ++i) {
		double val = 0.;
		myFile >> val;
		aVector.push_back(val);
	}

}

void fromFile(std::ifstream &myFile,  std::vector<Eigen::VectorXd> &aVector2) {

	size_t size = 0;
	myFile >> size;
	aVector2.clear();

	for(size_t i=0; i<size; ++i) {
		Eigen::VectorXd aVec;
		fromFile(myFile, aVec);
		aVector2.push_back(aVec);
	}

}

void fromFile(std::ifstream &myFile,  std::vector<Eigen::MatrixXd> &aTensor) {

	size_t N = 0;
	size_t rows = 0;
	size_t cols = 0;

	aTensor.clear();

	myFile >> N;

	if(N == 0) return;

	myFile >> rows;
	myFile >> cols;

	for(size_t k=0; k<N; ++k) {
		Eigen::MatrixXd aMatrix;
		aMatrix.resize(rows, cols);

		for(size_t i=0; i<rows; ++i) {
			for(size_t j=0; j<cols; ++j) {
				myFile >> aMatrix(i,j);
			}
		}

		aTensor.push_back(aMatrix);
	}

}

void fromFile(std::ifstream &myFile,  Eigen::VectorXd &aVector) {

	size_t size = 0;
	myFile >> size;
	aVector.resize(size);

	for(size_t i=0; i<size; ++i) {
		double val = 0.;
		myFile >> val;
		aVector(i) = val;
	}
}

void fromFile(std::ifstream &myFile,  Eigen::MatrixXd &aMatrix) {

	size_t rows = 0;
	size_t cols = 0;
	myFile >> rows;
	myFile >> cols;

	aMatrix.resize(rows, cols);

	for(size_t i=0; i<rows; ++i) {
		for(size_t j=0; j<cols; ++j) {
			myFile >> aMatrix(i,j);
		}
	}

}


void fromFile(std::ifstream &myFile, std::vector< Tensor::rbEventMap_t > &aVecSparseTensor) {

	aVecSparseTensor.clear();

	size_t nSparseTensor;
	myFile >> nSparseTensor;

	aVecSparseTensor.resize(nSparseTensor);

	for(size_t i=0; i<aVecSparseTensor.size(); ++i) {

		size_t nEntries;
		myFile >> nEntries;

		for(size_t j=0; j<nEntries; ++j) {
			std::vector<unsigned int> key(3, 0);
			myFile >> key[0];
			myFile >> key[1];
			myFile >> key[2];

			double val = 0.;
			myFile >> val;

			aVecSparseTensor[i][key] = val;
		}
	}

}


} /* namespace Parameters */
} /* namespace Serializer */
} /* namespace Utils */
