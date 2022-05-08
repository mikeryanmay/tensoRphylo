/*
 * IncLikelihoodApproximator.h
 *
 *  Created on: Sep 4, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_APPROXIMATOR_INCLIKELIHOODAPPROXIMATOR_H_
#define LIKELIHOOD_APPROXIMATOR_INCLIKELIHOODAPPROXIMATOR_H_

#include <boost/assign/list_of.hpp>
#include <vector>

#include "Factory.h"

#include "Specialized/SpecializedApproxManager.h"
#include "Specialized/BaseSpecialized.h"
#include "Specialized/BaseDense.h"
#include "Specialized/DenseUnobserved.h"
#include "Specialized/AdaptivePairedUS.h"

namespace Likelihood {
namespace Approximator {

#if defined(_OPENMP)

typedef enum {
	AUTO_TUNING=0,
	SEQUENTIAL_OPTIMIZED=1,
	SEQUENTIAL_BRANCHWISE=2,
	PARALLEL_OPTIMIZED=3,
	PARALLEL_BRANCHWISE=4
} approximatorVersion_t;

const std::vector<std::string> APPROXIMATOR_NAMES = boost::assign::list_of("APPROXIMATOR_AUTOTUNING")("APPROXIMATOR_SEQ_OPTIMIZED")("APPROXIMATOR_SEQ_BRANCHWISE")("APPROXIMATOR_PAR_OPTIMIZED")("APPROXIMATOR_PAR_BRANCHWISE");

#else

typedef enum {
	AUTO_TUNING=0,
	SEQUENTIAL_OPTIMIZED=1,
	SEQUENTIAL_BRANCHWISE=2
} approximatorVersion_t;

const std::vector<std::string> APPROXIMATOR_NAMES = boost::assign::list_of("APPROXIMATOR_AUTOTUNING")("APPROXIMATOR_SEQ_OPTIMIZED")("APPROXIMATOR_SEQ_BRANCHWISE");

#endif

}
}

#endif /* LIKELIHOOD_APPROXIMATOR_INCLIKELIHOODAPPROXIMATOR_H_ */
