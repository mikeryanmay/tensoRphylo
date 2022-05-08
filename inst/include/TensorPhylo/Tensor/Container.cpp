/*
 * Container.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#include "Container.h"

#include <iostream>
#include "IncTensor.h"

namespace Tensor {

Container::Container() {
}

Container::~Container() {
}

bool Container::isFullyDefined() const {
	return lambda != NULL &&
			    mu != NULL &&
			    phi != NULL &&
			    delta != NULL &&
			    omega != NULL &&
			    eta != NULL &&
				omega != NULL;
}


bool Container::areDimensionsMatching() const {
	bool vectorwise = lambda->getDimensions().size() == 1 &&
								   lambda->getDimensions()[0] > 0 &&
								   lambda->getDimensions() == mu->getDimensions() &&
								   lambda->getDimensions() == delta->getDimensions() &&
								   lambda->getDimensions() == phi->getDimensions();

	bool matrixwise = eta->getDimensions().size() == 2 &&
								   eta->getDimensions()[0] == eta->getDimensions()[1] &&
								   eta->getDimensions()[0] == lambda->getDimensions()[0];

	bool tensorwise = omega->getDimensions().size() == 3 &&
								   omega->getDimensions()[0] ==  omega->getDimensions()[1] &&
								   omega->getDimensions()[0] ==  omega->getDimensions()[2] &&
								   omega->getDimensions()[0] == lambda->getDimensions()[0];

	return vectorwise && matrixwise && tensorwise;
}

size_t Container::getNumberOfState() const {
	return lambda->getDimensions()[0];
}


void Container::setLambda(BaseTensorSharedPtr aPtrLambda) {
	lambda = aPtrLambda;
}

void Container::setMu(BaseTensorSharedPtr aPtrMu) {
	mu = aPtrMu;
}

void Container::setPhi(BaseTensorSharedPtr aPtrPhi) {
	phi = aPtrPhi;
}

void Container::setDelta(BaseTensorSharedPtr aPtrDelta) {
	delta = aPtrDelta;
}

void Container::setOmega(BaseTensorSharedPtr aPtrOmega) {
	omega = aPtrOmega;
}

void Container::setEta(BaseTensorSharedPtr aPtrEta) {
	eta = aPtrEta;
}

BaseTensorSharedPtr Container::getLambda() {
	return lambda;
}

BaseTensorSharedPtr Container::getMu() {
	return mu;
}

BaseTensorSharedPtr Container::getPhi() {
	return phi;
}

BaseTensorSharedPtr Container::getDelta() {
	return delta;
}

BaseTensorSharedPtr Container::getOmega() {
	return omega;
}

BaseTensorSharedPtr Container::getEta() {
	return eta;
}


Eigen::MatrixXd& Container::getEigenVecLambda( double t ) {
	assert(lambda);
	return lambda->getTensor(t).front();
}

Eigen::MatrixXd& Container::getEigenVecMu( double t ) {
	assert(mu);
	return mu->getTensor(t).front();
}

Eigen::MatrixXd& Container::getEigenVecPhi( double t ) {
	assert(phi);
	return phi->getTensor(t).front();
}

Eigen::MatrixXd& Container::getEigenVecDelta( double t ) {
	assert(delta);
	return delta->getTensor(t).front();
}

tensor_t& Container::getEigenTensorOmega( double t ) {
	return omega->getTensor(t);
}

Eigen::MatrixXd& Container::getEigenMatrixEta( double t ) {

	return eta->getTensor(t).front();
}


etaStructure_t Container::getEtaStructureType() const {

	EtaMatrix* etaMatrix = dynamic_cast<EtaMatrix*>(eta.get());
	assert(etaMatrix != NULL && "The tensor defined for eta is not of type Matrix.");

	return etaMatrix->getStructureType();
}


std::string Container::reportTensorInformations() const {
	std::stringstream ss;
	ss << "Tensors need to have less than " << 100.0*BaseTensor::NONZERO_COEFF_SPARSE_THRESHOLD << "% non-zero coefficients to be considered as sparse." << std::endl;
	ss << "Eta is: " << (eta->isSparse() ? "sparse; " : "dense; ") << (eta->isContantThroughTime() ? "homogeneous; " : "heterogeneous; ") << std::endl;
	ss << "Omega is: " << (omega->isSparse() ? "sparse; " : "dense; ") << (omega->isContantThroughTime() ? "homogeneous; " : "heterogeneous; ") << std::endl;
	return ss.str();
}

} /* namespace Tensor */
