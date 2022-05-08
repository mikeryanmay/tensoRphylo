/*
 * IncFwdStochasticMapping.h
 *
 *  Created on: May 11, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_INCFWDSTOCHASTICMAPPING_H_
#define LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_INCFWDSTOCHASTICMAPPING_H_
#include <boost/smart_ptr/shared_ptr.hpp>
#include <cstddef>
#include <map>
#include <vector>

namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

typedef size_t mapHistoriesKey_t;
typedef std::vector< std::pair<double, size_t> > mapHistoriesVal_t;
typedef std::map< mapHistoriesKey_t, mapHistoriesVal_t > mapHistories_t;

class SequentialStochasticMappingCPU;
typedef boost::shared_ptr<SequentialStochasticMappingCPU> StochasticMappingApproxSharedPtr;

typedef enum {REJECTION_SAMPLING_ALGO=0, DENSE_EULER_ALGO=1, DENSE_DOPRI_ALGO=2} stochasticMappingAlgo_t;


}
}
}

#endif /* LIKELIHOOD_APPROXIMATOR_STOCHASTICMAPPING_INCFWDSTOCHASTICMAPPING_H_ */
