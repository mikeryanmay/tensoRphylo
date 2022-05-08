/*
 * Factoy.cpp
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#include "Factory.h"

#include "Parameters/Container.h"
#include "Parameters/AsyncParameterContainer.h"
#include "Tensor/IncTensor.h"
#include "Parameters/IncParameterContainer.h"

namespace Tensor {

Factory::Factory() {
}

Factory::~Factory() {
}


Tensor::ContainerSharedPtr Factory::createContainerWithTimeHomogenousVectors(Parameters::ContainerSharedPtr ptrParameters) {

	BaseTensorSharedPtr ptrLambda( new Tensor::SimpleVector(ptrParameters->lambda) );
	BaseTensorSharedPtr ptrMu(     new Tensor::SimpleVector(ptrParameters->mu) );
	BaseTensorSharedPtr ptrPhi(    new Tensor::SimpleVector(ptrParameters->phi) );
	BaseTensorSharedPtr ptrDelta(  new Tensor::SimpleVector(ptrParameters->delta) );

	BaseTensorSharedPtr ptrEta(   new Tensor::SimpleMatrix(ptrParameters->eta) );

	BaseTensorSharedPtr ptrOmega( new Tensor::SimpleTensor(ptrParameters->omega) );

	Tensor::ContainerSharedPtr ptrContainer(new Tensor::Container);
	ptrContainer->setLambda(ptrLambda);
	ptrContainer->setMu(ptrMu);
	ptrContainer->setPhi(ptrPhi);
	ptrContainer->setDelta(ptrDelta);
	ptrContainer->setEta(ptrEta);
	ptrContainer->setOmega(ptrOmega);

	return ptrContainer;
}



Tensor::ContainerSharedPtr Factory::createContainerWithTimeHeterogenousVectors(Parameters::AsyncContainerSharedPtr ptrParameters) {

	assert(!ptrParameters->vecLambda.empty() || !ptrParameters->vecMu.empty());
	size_t nState = 0;
	if(!ptrParameters->vecLambda.empty()) {
		nState = ptrParameters->vecLambda.front().size();
	} else {
		nState = ptrParameters->vecMu.front().size();
	}
	assert(nState > 0);

	if(ptrParameters->vecLambda.empty()) {
		ptrParameters->vecLambda.push_back(Eigen::VectorXd::Zero(nState));
	}
	BaseTensorSharedPtr ptrLambda( new Tensor::TimeHeterogenousVector(ptrParameters->timesLambda, ptrParameters->vecLambda) );

	if(ptrParameters->vecMu.empty()) {
		ptrParameters->vecMu.push_back(Eigen::VectorXd::Zero(nState));
	}
	BaseTensorSharedPtr ptrMu(     new Tensor::TimeHeterogenousVector(ptrParameters->timesMu, ptrParameters->vecMu) );

	if(ptrParameters->vecPhi.empty()) {
		ptrParameters->vecPhi.push_back(Eigen::VectorXd::Zero(nState));
	}
	BaseTensorSharedPtr ptrPhi(    new Tensor::TimeHeterogenousVector(ptrParameters->timesPhi, ptrParameters->vecPhi) );

	if(ptrParameters->vecDelta.empty()) {
		ptrParameters->vecDelta.push_back(Eigen::VectorXd::Zero(nState));
	}
	BaseTensorSharedPtr ptrDelta(  new Tensor::TimeHeterogenousVector(ptrParameters->timesDelta, ptrParameters->vecDelta) );

	if(ptrParameters->vecEta.empty()) {
		ptrParameters->vecEta.push_back(Eigen::MatrixXd::Zero(nState, nState));
	}
	BaseTensorSharedPtr ptrEta(   new Tensor::TimeHeterogenousMatrix(ptrParameters->timesEta, ptrParameters->vecEta) );

	BaseTensorSharedPtr ptrOmega( new Tensor::TimeHeterogenousTensor(nState, ptrParameters->timesOmega, ptrParameters->vecOmega) );

	Tensor::ContainerSharedPtr ptrContainer(new Tensor::Container);
	ptrContainer->setLambda(ptrLambda);
	ptrContainer->setMu(ptrMu);
	ptrContainer->setPhi(ptrPhi);
	ptrContainer->setDelta(ptrDelta);
	ptrContainer->setEta(ptrEta);
	ptrContainer->setOmega(ptrOmega);

	return ptrContainer;
}

} /* namespace Tensor */
