/*
 * Container.cpp
 *
 *  Created on: Sep 3, 2019
 *      Author: xaviermeyer
 */

#include "Container.h"

#include "Definition.h"
#include "Parameters/IncParameterContainer.h"

namespace SynchronousEvents {

Container::Container(Parameters::ContainerSharedPtr aPtrParameters) {

	ptrMassSpeciation.reset(new Definition(aPtrParameters->massSpeciationTimes, aPtrParameters->massSpeciationProb, emptyVector));
	ptrMassExtinction.reset(new Definition(aPtrParameters->massExtinctionTimes, aPtrParameters->massExtinctionProb, aPtrParameters->massExtinctionStateChangeProb));
	ptrMassSampling.reset(new Definition(aPtrParameters->massSamplingTimes, aPtrParameters->massSamplingProb, emptyVector));
	ptrMassDestrSampling.reset(new Definition(aPtrParameters->massDestrSamplingTimes, aPtrParameters->massDestrSamplingProb, emptyVector));

}

Container::Container(Parameters::SyncContainerSharedPtr aPtrParameters) {

	ptrMassSpeciation.reset(new Definition(aPtrParameters->massSpeciationTimes, aPtrParameters->massSpeciationProb, emptyVector));
	ptrMassExtinction.reset(new Definition(aPtrParameters->massExtinctionTimes, aPtrParameters->massExtinctionProb, aPtrParameters->massExtinctionStateChangeProb));
	ptrMassSampling.reset(new Definition(aPtrParameters->massSamplingTimes, aPtrParameters->massSamplingProb, emptyVector));
	ptrMassDestrSampling.reset(new Definition(aPtrParameters->massDestrSamplingTimes, aPtrParameters->massDestrSamplingProb, emptyVector));

}

Container::~Container() {
}

DefinitionSharedPtr Container::getPtrMassSpeciation() {
	return ptrMassSpeciation;
}

DefinitionSharedPtr Container::getPtrMassExtinction() {
	return ptrMassExtinction;
}

DefinitionSharedPtr Container::getPtrMassSampling() {
	return ptrMassSampling;
}

DefinitionSharedPtr Container::getPtrMassDestrSampling() {
	return ptrMassDestrSampling;
}



} /* namespace SynchronousEvents */
