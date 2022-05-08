/*
 * SegmentStatistics.cpp
 *
 *  Created on: May 11, 2020
 *      Author: meyerx
 */

#include "SegmentStatistics.h"

#include <iostream>


namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

SegmentStatistics::SegmentStatistics(const size_t aNStates) : N_STATE(aNStates) {

	count = 0;
	nTransitions = 0.;
	timeInStates = Eigen::VectorXd::Zero(N_STATE);
	cntTransitions = Eigen::MatrixXd::Zero(N_STATE, N_STATE);

	segmentLength = 0.;
}

SegmentStatistics::~SegmentStatistics() {
}


void SegmentStatistics::registerSegmentHistory(const mapHistoriesVal_t &branchHistory) {

	assert(branchHistory.size() >= 1);
	nTransitions += branchHistory.size()-1;

	for(size_t iE=0; iE<branchHistory.size(); ++iE) {
		size_t curState = branchHistory[iE].second;
		double curDuration = branchHistory[iE].first;
		timeInStates(curState) += fabs(curDuration);

		if(iE+1 < branchHistory.size()) {
			cntTransitions(curState, branchHistory[iE+1].second) += 1.;
		}

		if(count == 0) {
			segmentLength += curDuration;
		}
	}
	count++;
}

size_t SegmentStatistics::getCount() const {
	return count;
}

double SegmentStatistics::getSegmentLength() const {
	return segmentLength;
}

double SegmentStatistics::getAverageNTransitions() const {
	return nTransitions / count;
}

Eigen::VectorXd SegmentStatistics::getStateFrequency() const {
	Eigen::VectorXd stateFreq = timeInStates;
	stateFreq /= stateFreq.sum();
	return stateFreq;
}

Eigen::MatrixXd SegmentStatistics::getTransitionFrequency() const {
	Eigen::MatrixXd transitionFrequencies = cntTransitions;
	transitionFrequencies /= transitionFrequencies.sum();
	return transitionFrequencies;
}


std::string SegmentStatistics::summarizedStatsToString() const {
	std::stringstream ss;
	ss << "Mean nb transitions : " << std::endl;
	ss << getAverageNTransitions() << std::endl;
	ss << "State frequency : " << std::endl;
	ss << getStateFrequency() << std::endl;
	ss << "Transitions : " << std::endl;
	ss << getTransitionFrequency() << std::endl;

	return ss.str();
}


bool areSegmentStatisticsSimilar(const SegmentStatistics &segmentA, const SegmentStatistics &segmentB, double tolerance) {
	/*std::cout << "Segment A : " << std::endl << segmentA.summarizedStatsToString() << std::endl;
	std::cout << "Segment B : " << std::endl << segmentB.summarizedStatsToString() << std::endl;
	std::cout << "-----------------------------------------" << std::endl;*/

	if(segmentA.getSegmentLength() == 0. && segmentB.getSegmentLength() == 0.) return true;

	bool areSimilar = true;
	areSimilar = areSimilar && fabs(segmentA.getAverageNTransitions() - segmentB.getAverageNTransitions()) < tolerance;
	areSimilar = areSimilar && (segmentA.getStateFrequency() - segmentB.getStateFrequency()).norm()/segmentA.getStateFrequency().size() < tolerance;
	if(segmentA.getAverageNTransitions()*segmentA.getCount() >= 10. && segmentB.getAverageNTransitions()*segmentB.getCount() >= 10.) {
		areSimilar = areSimilar && (segmentA.getTransitionFrequency() - segmentB.getTransitionFrequency()).norm()/segmentA.getTransitionFrequency().size() < tolerance;
	}

	if(!areSimilar) {
		std::cout << fabs(segmentA.getAverageNTransitions() - segmentB.getAverageNTransitions()) << std::endl;
		std::cout << (segmentA.getStateFrequency() - segmentB.getStateFrequency()).norm() << std::endl;
		std::cout << (segmentA.getTransitionFrequency() - segmentB.getTransitionFrequency()).norm() << std::endl;
		std::cout << "Segment A : " << std::endl << segmentA.summarizedStatsToString() << std::endl;
		std::cout << "Segment B : " << std::endl << segmentB.summarizedStatsToString() << std::endl;
		std::cout << "-----------------------------------------" << std::endl;
	}

	return areSimilar;
}

} /* namespace StochasticMapping */
} /* namespace Approximator */
} /* namespace Likelihood */
