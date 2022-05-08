/*
 * IncLikelihoodApproximator.h
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_INCFWDLIKELIHOODAPPROXIMATOR_H_
#define LIKELIHOOD_APPROXIMATOR_INCFWDLIKELIHOODAPPROXIMATOR_H_

#include <boost/smart_ptr/shared_ptr.hpp>
#include "StochasticMapping/IncFwdStochasticMapping.h"

namespace Likelihood {
namespace Approximator {

class BaseApproximator;
typedef boost::shared_ptr<BaseApproximator> ApproximatorSharedPtr;

namespace Specialized {

class BaseSpecialized;
typedef boost::shared_ptr<BaseSpecialized> BaseSpecializedSharedPtr;

class BaseDense;
typedef boost::shared_ptr<BaseDense> BaseDenseSharedPtr;

class DenseUnobserved;
typedef boost::shared_ptr<DenseUnobserved> DenseUnobservedSharedPtr;

class AdaptivePairedUS;
typedef boost::shared_ptr<AdaptivePairedUS> AdaptivePairedUSSharedPtr;

class SpecializedApproxManager;
typedef boost::shared_ptr<SpecializedApproxManager> SpecializedApproxManagerSharedPtr;

typedef enum {EXTANT_PROCESS, FULL_PROCESS} processType_t;

}

}
}

#endif /* LIKELIHOOD_APPROXIMATOR_INCFWDLIKELIHOODAPPROXIMATOR_H_ */
