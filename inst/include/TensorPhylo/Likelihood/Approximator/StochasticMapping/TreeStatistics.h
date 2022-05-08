/*
 * TreeStatistics.h
 *
 *  Created on: May 11, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_TREESTATISTICS_H_
#define LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_TREESTATISTICS_H_

#include <Eigen/Core>
#include "Interface/DistributionHandler.h"
#include "IncFwdStochasticMapping.h"
#include "SegmentStatistics.h"

namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

class TreeStatistics {
public:
	TreeStatistics(size_t aNStates, const TensorPhylo::Interface::vecHistories_t &vecHistories);
	~TreeStatistics();

	SegmentStatistics& getSegmentStatistics(size_t iB);

	std::pair<double, double> getObservedTransitionFrequenciesStats() const;
	const Eigen::MatrixXd& getObservedTransitionMatrix() const;

private:

	const size_t N_STATES;

	typedef std::map<size_t, SegmentStatistics> mapSegmentStats_t;
	mapSegmentStats_t mapSegmentStats;
	Eigen::MatrixXd observedTransisionMatrix;


	void registerHistory(const mapHistories_t &history);

	friend bool areTreeStatisticsSimilar(const TreeStatistics &treeStatsA, const TreeStatistics &treeStatsB, double tolerance);

};

bool areTreeStatisticsSimilar(const TreeStatistics &treeStatsA, const TreeStatistics &treeStatsB, double tolerance);

} /* namespace StochasticMapping */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_TREESTATISTICS_H_ */
