/*
 * SegmentStatistics.h
 *
 *  Created on: May 11, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_SEGMENTSTATISTICS_H_
#define LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_SEGMENTSTATISTICS_H_

#include <Eigen/Core>
#include "IncFwdStochasticMapping.h"

namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

class SegmentStatistics {
public:
	SegmentStatistics(const size_t aNStates);
	~SegmentStatistics();

	void registerSegmentHistory(const mapHistoriesVal_t &segmentHistory);

	size_t getCount() const;
	double getSegmentLength() const;
	double getAverageNTransitions() const;
	Eigen::VectorXd getStateFrequency() const;
	Eigen::MatrixXd getTransitionFrequency() const;

	std::string summarizedStatsToString() const;

private:

	const size_t N_STATE;

	double segmentLength;
	size_t count;

	double nTransitions;
	Eigen::VectorXd timeInStates;
	Eigen::MatrixXd cntTransitions;

};

bool areSegmentStatisticsSimilar(const SegmentStatistics &segmentA, const SegmentStatistics &segmentB, double tolerance);

} /* namespace StochasticMapping */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_SEGMENTSTATISTICS_H_ */
