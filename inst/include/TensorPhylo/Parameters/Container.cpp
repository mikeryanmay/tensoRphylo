/*
 * ParametersContainer.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#include "Container.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>

namespace Parameters {

Container::Container() {

	applyTreeLikCorrection = false;
	intLikApproximator = 1;
	nThreads = 1;

	intScheme = 0;
	condType = 0;
	deltaT = 0.001;

}

Container::~Container() {
}

void Container::readFromFile(const std::string &filename) {
	std::ifstream myFile;
	myFile.open(filename.c_str(), std::ifstream::in);
	assert(myFile.good());

	std::string line;
	while(std::getline(myFile, line)) {

		boost::algorithm::trim(line);
		if(line.empty()) continue;

		if(line =="applyTreeLikCorrection") {
			std::string strApplyTLC;
			myFile >> strApplyTLC;
			boost::algorithm::to_lower<std::string>(strApplyTLC);
			assert((strApplyTLC == "true" || strApplyTLC == "false") && "The value for 'likTreeCorrection' must be 'true' or 'false'.");
			applyTreeLikCorrection = strApplyTLC == "true";
		} else if(line =="intLikApproximator") {
			myFile >> intLikApproximator;
		} else if(line == "nThreads") {
			myFile >> nThreads;
		} else if(line =="intScheme") {
			myFile >> intScheme;
		} else if(line =="condType") {
			myFile >> condType;
		} else if(line =="deltaT") {
			myFile >> deltaT;
		} else if(line =="rootPrior") {
			fromFile(myFile, rootPrior);
		} else if(line =="lambda") {
			fromFile(myFile, lambda);
		} else if(line =="mu") {
			fromFile(myFile, mu);
		} else if(line =="phi") {
			fromFile(myFile, phi);
		} else if(line =="delta") {
			fromFile(myFile, delta);
		} else if(line =="eta") {
			fromFile(myFile, eta);
		} else if(line =="omega") {
			fromFile(myFile, omega);
		} else if(line =="massSpeciationTimes") {
			fromFile(myFile, massSpeciationTimes);
		} else if(line =="massSpeciationProb") {
			fromFile(myFile, massSpeciationProb);
		} else if(line =="massExtinctionTimes") {
			fromFile(myFile, massExtinctionTimes);
		} else if(line =="massExtinctionProb") {
			fromFile(myFile, massExtinctionProb);
		} else if(line =="massExtinctionStateChangeProb") {
			fromFile(myFile, massExtinctionStateChangeProb);
		} else if(line =="massSamplingTimes") {
			fromFile(myFile, massSamplingTimes);
		} else if(line =="massSamplingProb") {
			fromFile(myFile, massSamplingProb);
		} else if(line =="massDestrSamplingTimes") {
			fromFile(myFile, massDestrSamplingTimes);
		} else if(line =="massDestrSamplingProb") {
			fromFile(myFile, massDestrSamplingProb);
		} else if(line =="synchMonitoring") {
			fromFile(myFile, synchMonitoring);
		} else {
			std::cerr << "Parameters reader error at string: " << line << std::endl;
		}

	}

	myFile.close();


}

void Container::writeToFile(const std::string &filename) {

	std::ofstream myFile;
	myFile.open(filename.c_str(), std::ofstream::out);
	assert(myFile.good());

	myFile << "applyTreeLikCorrection" << std::endl;
	if(applyTreeLikCorrection) {
		myFile << "true" << std::endl << std::endl;
	} else {
		myFile << "false" << std::endl << std::endl;
	}

	myFile << "intLikApproximator" << std::endl;
	myFile << intLikApproximator << std::endl << std::endl;

	myFile << "nThreads" << std::endl;
	myFile << nThreads << std::endl << std::endl;

	myFile << "intScheme" << std::endl;
	myFile << intScheme << std::endl << std::endl;

	myFile << "condType" << std::endl;
	myFile << condType << std::endl << std::endl;

	myFile << "deltaT" << std::endl;
	myFile << deltaT << std::endl << std::endl;

	myFile << "rootPrior" << std::endl;
	myFile << toString(rootPrior) << std::endl ;

	myFile << "lambda" << std::endl;
	myFile << toString(lambda) << std::endl;

	myFile << "mu" << std::endl;
	myFile << toString(mu) << std::endl;

	myFile << "phi" << std::endl;
	myFile << toString(phi) << std::endl;

	myFile << "delta" << std::endl;
	myFile << toString(delta) << std::endl;

	myFile << "eta" << std::endl;
	myFile << toString(eta) << std::endl;

	myFile << "omega" << std::endl;
	myFile << toString(omega) << std::endl;

	myFile << "massSpeciationTimes" << std::endl;
	myFile << toString(massSpeciationTimes) << std::endl;

	myFile << "massSpeciationProb" << std::endl;
	myFile << toString(massSpeciationProb) << std::endl;

	myFile << "massExtinctionTimes" << std::endl;
	myFile << toString(massExtinctionTimes) << std::endl;

	myFile << "massExtinctionProb" << std::endl;
	myFile << toString(massExtinctionProb) << std::endl;

	myFile << "massExtinctionStateChangeProb" << std::endl;
	myFile << toString(massExtinctionStateChangeProb) << std::endl;

	myFile << "massSamplingTimes" << std::endl;
	myFile << toString(massSamplingTimes) << std::endl;

	myFile << "massSamplingProb" << std::endl;
	myFile << toString(massSamplingProb) << std::endl;

	myFile << "massDestrSamplingTimes" << std::endl;
	myFile << toString(massDestrSamplingTimes) << std::endl;

	myFile << "massDestrSamplingProb" << std::endl;
	myFile << toString(massDestrSamplingProb) << std::endl;

	myFile << "synchMonitoring" << std::endl;
	myFile << toString(synchMonitoring) << std::endl;

	myFile.close();

}


std::string Container::toString(const std::vector<double> &aVector) {
	std::stringstream ss;
	ss << aVector.size() << std::endl;
	for(size_t i=0; i<aVector.size(); ++i) {
		ss << aVector[i];
		if(i != aVector.size()-1) ss << " ";
	}
	if(!aVector.empty()) ss << std::endl;

	return ss.str();
}

std::string Container::toString(const std::vector<Eigen::VectorXd> &aVector2) {
	std::stringstream ss;
	ss << aVector2.size() << std::endl;
	for(size_t i=0; i<aVector2.size(); ++i) {
		ss << toString(aVector2[i]);
	}
	return ss.str();
}

std::string Container::toString(const std::vector<Eigen::MatrixXd> &aTensor) {
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

std::string Container::toString(const Eigen::VectorXd &aVector) {
	std::stringstream ss;
	ss << aVector.size() << std::endl;
	for(size_t i=0; i<(size_t)aVector.size(); ++i) {
		ss << aVector(i);
		if(i != (size_t)aVector.size()-1) ss << " ";
	}
	ss << std::endl;

	return ss.str();
}

std::string Container::toString(const Eigen::MatrixXd &aMatrix) {
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

void Container::fromFile(std::ifstream &myFile,  std::vector<double> &aVector) {

	size_t size = 0;
	myFile >> size;
	aVector.clear();

	for(size_t i=0; i<size; ++i) {
		double val = 0.;
		myFile >> val;
		aVector.push_back(val);
	}

}

void Container::fromFile(std::ifstream &myFile,  std::vector<Eigen::VectorXd> &aVector2) {

	size_t size = 0;
	myFile >> size;
	aVector2.clear();

	for(size_t i=0; i<size; ++i) {
		Eigen::VectorXd aVec;
		fromFile(myFile, aVec);
		aVector2.push_back(aVec);
	}

}

void Container::fromFile(std::ifstream &myFile,  std::vector<Eigen::MatrixXd> &aTensor) {

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

void Container::fromFile(std::ifstream &myFile,  Eigen::VectorXd &aVector) {

	size_t size = 0;
	myFile >> size;
	aVector.resize(size);

	for(size_t i=0; i<size; ++i) {
		double val = 0.;
		myFile >> val;
		aVector(i) = val;
	}
}

void Container::fromFile(std::ifstream &myFile,  Eigen::MatrixXd &aMatrix) {

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


} /* namespace Parameters */
