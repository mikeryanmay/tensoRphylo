/*
 * TreeStatistics.cpp
 *
 *  Created on: May 11, 2020
 *      Author: meyerx
 */

#include "TreeStatistics.h"

#include <iostream>


namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

TreeStatistics::TreeStatistics(size_t aNStates, const TensorPhylo::Interface::vecHistories_t &vecHistories) :
		N_STATES(aNStates) {

	for(size_t iH=0; iH<vecHistories.size(); ++iH) {
		registerHistory(vecHistories[iH]);
	}


	bool first = true;
	for(mapSegmentStats_t::iterator itM = mapSegmentStats.begin(); itM != mapSegmentStats.end(); ++itM) {

		if(itM->second.getAverageNTransitions()*itM->second.getCount() >= 10.) {
			if(first) {
				observedTransisionMatrix = itM->second.getTransitionFrequency();
				first = false;
			} else {
				observedTransisionMatrix += itM->second.getTransitionFrequency();
			}
		}
	}
	observedTransisionMatrix /= observedTransisionMatrix.sum();
}

TreeStatistics::~TreeStatistics() {
}

SegmentStatistics& TreeStatistics::getSegmentStatistics(size_t iB) {
	mapSegmentStats_t::iterator itFind = mapSegmentStats.find(iB);
	assert(itFind != mapSegmentStats.end());
	return itFind->second;
}

const Eigen::MatrixXd& TreeStatistics::getObservedTransitionMatrix() const {
	return observedTransisionMatrix;
}

std::pair<double, double> TreeStatistics::getObservedTransitionFrequenciesStats() const {

	std::pair<double, double> transitionFreqFirstMoments;

	size_t cnt = 0;
	Eigen::VectorXd values(N_STATES*(N_STATES-1));
	for(size_t i=0; i<(size_t)observedTransisionMatrix.rows(); ++i) {
		for(size_t j=0; j<(size_t)observedTransisionMatrix.cols(); ++j) {
			if(i==j) continue;
			values(cnt) = observedTransisionMatrix(i,j);
			cnt++;
		}
	}

	transitionFreqFirstMoments.first = values.mean();
	for(size_t i=0; i<(size_t)values.size(); ++i) {
		transitionFreqFirstMoments.second = std::pow(values(i) - transitionFreqFirstMoments.first, 2);
	}
	transitionFreqFirstMoments.second /= (values.size()-1);
	transitionFreqFirstMoments.second = std::sqrt(transitionFreqFirstMoments.second);

	return transitionFreqFirstMoments;
}

void TreeStatistics::registerHistory(const mapHistories_t &history) {
	bool first = mapSegmentStats.empty();
	for(mapHistories_t::const_iterator itH = history.begin(); itH != history.end(); ++itH) {
		if(first) {
			SegmentStatistics segStats(N_STATES);
			segStats.registerSegmentHistory(itH->second);
			mapSegmentStats.insert(std::make_pair(itH->first, segStats));
		} else {
			mapSegmentStats_t::iterator itFind = mapSegmentStats.find(itH->first);
			assert(itFind != mapSegmentStats.end());
			itFind->second.registerSegmentHistory(itH->second);
		}
	}
}

bool areTreeStatisticsSimilar(const TreeStatistics &treeStatsA, const TreeStatistics &treeStatsB, double tolerance) {

	bool areSimilar = true;

	for(TreeStatistics::mapSegmentStats_t::const_iterator it = treeStatsA.mapSegmentStats.begin(); it != treeStatsA.mapSegmentStats.end(); ++it) {
		TreeStatistics::mapSegmentStats_t::const_iterator itFind = treeStatsB.mapSegmentStats.find(it->first);
		if(itFind == treeStatsB.mapSegmentStats.end()) {
			std::cout << "Missing branch : " << it->first << std::endl;
			return false;
		}

		areSimilar = areSimilar && areSegmentStatisticsSimilar(it->second, itFind->second, tolerance);
	}

	return areSimilar;

}

} /* namespace StochasticMapping */
} /* namespace Approximator */
} /* namespace Likelihood */
