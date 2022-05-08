/*
 * BaseStochasticMappingApproximator.h
 *
 *  Created on: May 8, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_BASESTOCHASTICMAPPINGAPPROXIMATOR_H_
#define LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_BASESTOCHASTICMAPPINGAPPROXIMATOR_H_

#include "IncFwdStochasticMapping.h"

namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

class BaseStochasticMappingApproximator {
public:
	BaseStochasticMappingApproximator();
	virtual ~BaseStochasticMappingApproximator();

	virtual mapHistories_t drawHistory() = 0;
	virtual mapHistories_t drawHistoryAndComputeRates(std::vector<double>& averageLambda, std::vector<double>& averageMu, std::vector<double>& averagePhi, std::vector<double>& averageDelta, std::vector<long>& numChanges) = 0;
	virtual mapHistories_t drawAncestralStates() = 0;

};

} /* namespace StochasticMapping */
} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_BASESTOCHASTICMAPPINGAPPROXIMATOR_H_ */
