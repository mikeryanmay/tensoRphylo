/*
 * BaseApproximator.cpp
 *
 *  Created on: Aug 25, 2019
 *      Author: meyerx
 */

#include "BaseApproximator.h"

#include "../../Data/Structure/Node.h"
#include "Data/Structure/Tree.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Likelihood/ConditionTypes/ConditionType.h"
#include "Likelihood/Kernels/CPU/Misc/QuasistationaryFrequency.h"

namespace Likelihood {
namespace Approximator {

const double BaseApproximator::DEFAULT_DELTA_T = 5.E-1;
const double BaseApproximator::DEFAULT_ABS_TOLERANCE = 1.E-7;
const double BaseApproximator::DEFAULT_REL_TOLERANCE = 1.E-7;

BaseApproximator::BaseApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
								   Conditions::conditionalProbability_t aConditionType,
									 Phylogeny::Data::ContainerSharedPtr aPtrData,
								   Scheduler::SchedulerSharedPtr aPtrScheduler,
								   SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
								   Tensor::ContainerSharedPtr aPtrTensorCont) :
											intScheme(aIntScheme),
											conditionType(aConditionType),
											ptrData(aPtrData), ptrScheduler(aPtrScheduler),
											ptrSyncEventsCont(aPtrSyncEventsCont), ptrTensorCont(aPtrTensorCont) {


	deltaT = DEFAULT_DELTA_T;
	applyTreeCorrection = true;
	condCompatibilityMode = false;
	useQuasistationaryFrequency = false;
	logLikelihood = 0.;
}

BaseApproximator::~BaseApproximator() {
}

void BaseApproximator::enableTreeLikelihoodCorrection() {
	applyTreeCorrection = true;
}

void BaseApproximator::disableTreeLikelihoodCorrection() {
	applyTreeCorrection = false;
}

void BaseApproximator::setConditionalProbabilityCompatibilityMode(bool setActive) {
	condCompatibilityMode = setActive;
}

void BaseApproximator::setQuasistationaryFrequencyMode(bool setActive) {
	useQuasistationaryFrequency = setActive;
}

Eigen::VectorXd BaseApproximator::getRootFrequency(double t) {
	Eigen::VectorXd sf;
	if ( useQuasistationaryFrequency ) {
		sf = Likelihood::Kernels::CPU::computeQuasiStationaryFrequency(ptrTensorCont, t);
	} else {
		sf = priorStateProbability;
	}
	return sf;
}

void BaseApproximator::setPriorStateProbability(const Eigen::VectorXd &aPriorStateProbability) {
	priorStateProbability = aPriorStateProbability;
}

double BaseApproximator::approximateLikelihood() {
	return exp(approximateLogLikelihood());
}

void BaseApproximator::orderProbes() {

	// Consolidate probes
	if(!vecProbesState.empty()) {
		// Order by time and merge similar times
		std::vector<Likelihood::Monitor::ProbeState> tmpVecPS = vecProbesState;
		vecProbesState.clear();

		std::sort(tmpVecPS.begin(), tmpVecPS.end());
		for(size_t iP=0; iP<tmpVecPS.size(); ++iP) {
			if(vecProbesState.empty() || vecProbesState.back().time != tmpVecPS[iP].time) {
				vecProbesState.push_back(tmpVecPS[iP]);
			} else {
				vecProbesState.back().vecP.push_back(tmpVecPS[iP].vecP.front());
				vecProbesState.back().vecIdEdge.push_back(tmpVecPS[iP].vecIdEdge.front());
				vecProbesState.back().scalingFactor.push_back(tmpVecPS[iP].scalingFactor.front());
			}
		}

		// Order by edges
		for(size_t iP=0; iP<vecProbesState.size(); ++iP) {

			// populate with arg
			std::vector< std::pair<long int, size_t> > sortedEdges(vecProbesState[iP].vecIdEdge.size());
			for(size_t iE=0; iE<vecProbesState[iP].vecIdEdge.size(); ++iE) {
				sortedEdges[iE] = std::make_pair(vecProbesState[iP].vecIdEdge[iE], iE);
			}

			// edge + argsort
			std::sort(sortedEdges.begin(), sortedEdges.end());

			// sort
			Likelihood::Monitor::ProbeState tmpProbe = vecProbesState[iP];
			for(size_t iE=0; iE<vecProbesState[iP].vecIdEdge.size(); ++iE) {
				vecProbesState[iP].vecIdEdge[iE] = tmpProbe.vecIdEdge[sortedEdges[iE].second];
				vecProbesState[iP].scalingFactor[iE] = tmpProbe.scalingFactor[sortedEdges[iE].second];
				vecProbesState[iP].vecP[iE] = tmpProbe.vecP[sortedEdges[iE].second];
			}
		}
	}
}


const std::vector<Likelihood::Monitor::ProbeState>& BaseApproximator::getObservedProbesState() const {
	return vecProbesState;
}


const std::vector<double>& BaseApproximator::getIntegrationTimes() const{
	return integrationTimes;
}

bool BaseApproximator::areEventsPossible() {

	// Check the event validity
	for(size_t iE=0; iE<ptrScheduler->getEvents().size(); ++iE) {
		if(!ptrScheduler->getEvents()[iE]->isEventPossible()) { // If this event is impossible
			return false; // We are done
		}
	}

	// All events are possible
	return true;
}


} /* namespace Approximator */
} /* namespace Likelihood */
