/*
 * Definition.cpp
 *
 *  Created on: Sep 3, 2019
 *      Author: xaviermeyer
 */

#include "Definition.h"

namespace SynchronousEvents {

Definition::Definition(std::vector<double> &aEventTimes, std::vector<Eigen::VectorXd> &aEventProbabilities, std::vector<Eigen::MatrixXd> &aStateChangeProb) :
	eventTimes(aEventTimes), eventProbabilities(aEventProbabilities), stateChangeProb(aStateChangeProb) {

}

Definition::~Definition() {
}

bool Definition::hasEvents() const {
	return eventTimes.size() > 0;
}


bool Definition::hasEvent(double t) const {
	for(size_t i=0; i<eventTimes.size(); ++i) {
		bool areSameUpToEpsilon = (std::fabs(t -eventTimes[i]) <= 2.*std::numeric_limits<double>::epsilon());
		if(areSameUpToEpsilon) {
			return true;
		}
	}
	return false;
}


bool Definition::areStateChangeInvolved() const {
	return !stateChangeProb.empty();
}

std::vector<double>& Definition::getEventTimes() {
	return eventTimes;
}

Eigen::VectorXd& Definition::getEventProbability(double t) {
	//std::vector<double>::iterator itFind = std::find(eventTimes.begin(), eventTimes.end(), t);
	//assert(itFind != eventTimes.end());
	size_t idx = 0;
	bool found = false;
	for(size_t i=0; i<eventTimes.size(); ++i) {
		bool areSameUpToEpsilon = (std::fabs(t -eventTimes[i]) <= 2.*std::numeric_limits<double>::epsilon());
		if(areSameUpToEpsilon) {
			found = true;
			idx = i;
			break;
		}
	}
	assert(found);

	//size_t idx = std::distance(eventTimes.begin(), itFind);
	return eventProbabilities[idx];
}

Eigen::MatrixXd& Definition::getStateChangeProbability(double t) {
	assert(areStateChangeInvolved());
	//std::vector<double>::iterator itFind = std::find(eventTimes.begin(), eventTimes.end(), t);
	//assert(itFind != eventTimes.end());
	size_t idx = 0;
	bool found = false;
	for(size_t i=0; i<eventTimes.size(); ++i) {
		bool areSameUpToEpsilon = (std::fabs(t -eventTimes[i]) <= 2.*std::numeric_limits<double>::epsilon());
		if(areSameUpToEpsilon) {
			found = true;
			idx = i;
			break;
		}
	}
	assert(found);

	//size_t idx = std::distance(eventTimes.begin(), itFind);
	return stateChangeProb[idx];
}

} /* namespace SynchronousEvent */
