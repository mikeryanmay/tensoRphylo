/*
 * Container.h
 *
 *  Created on: Sep 3, 2019
 *      Author: xaviermeyer
 */

#ifndef SYNCHRONOUSEVENT_CONTAINER_H_
#define SYNCHRONOUSEVENT_CONTAINER_H_

#include "IncFwdSynchronousEvents.h"
#include "Parameters/IncFwdParameterContainer.h"

#include <vector>
#include <Eigen/Core>

namespace SynchronousEvents {

class Container {
public:
	Container(Parameters::ContainerSharedPtr aPtrParameters);
	Container(Parameters::SyncContainerSharedPtr aPtrParameters);
	~Container();

	DefinitionSharedPtr getPtrMassSpeciation();
	DefinitionSharedPtr getPtrMassExtinction();
	DefinitionSharedPtr getPtrMassSampling();
	DefinitionSharedPtr getPtrMassDestrSampling();

private:

	std::vector<Eigen::MatrixXd> emptyVector;
	DefinitionSharedPtr ptrMassSpeciation, ptrMassExtinction, ptrMassSampling, ptrMassDestrSampling;

};

} /* namespace SynchronousEvents */

#endif /* SYNCHRONOUSEVENT_CONTAINER_H_ */
