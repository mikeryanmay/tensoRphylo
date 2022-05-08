/*
 * ConditionType.h
 *
 *  Created on: Sep 11, 2019
 *      Author: mrmay
 */

#ifndef LIKELIHOOD_CONDITIONTYPES_CONDITIONTYPE_H_
#define LIKELIHOOD_CONDITIONTYPES_CONDITIONTYPE_H_

#include <boost/assign/list_of.hpp>
#include <vector>
#include <cstddef>

namespace Likelihood {
namespace Conditions {

typedef enum {
	TIME=0,
	ROOT_SURVIVAL=1, 			// Require unobserved - no sampling 			(u_hat)
	ROOT_MRCA=2, 				// Require unobserved - no sampling 			(u_hat)
	ROOT_SAMPLING=3, 			// Require unobserved 							(u)
	STEM_SURVIVAL=4, 			// Require unobserved - no sampling 			(u_hat)
	STEM_ONE_SAMPLE=5, 			// Require unobserved 							(u)
	STEM_TWO_EXT_SAMPLES=6, 	// Require unobserved, singleton - no sampling 	(u_hat, s_hat)
	STEM_TWO_SAMPLES=7, 		// Require unobserved, singleton				(u, s)
	ROOT_SAMPLING_AND_MRCA=8 	// Require unobserved, singleton				(u)
} conditionalProbability_t;

const std::vector<std::string> CONDITION_NAMES = boost::assign::list_of("CONDITION_TIME")("CONDITION_ROOT_SURVIVAL")("CONDITION_ROOT_MRCA")("CONDITION_ROOT_SAMPLING")("CONDITION_STEM_SURVIVAL")("CONDITION_STEM_ONE_SAMPLE")("CONDITION_STEM_TWO_EXTANT")("CONDITION_STEM_TWO_SAMPLES")("CONDITION_ROOT_SAMPLING_AND_MRCA");

conditionalProbability_t intToConditionalProbabilityType(int myInt);

size_t getNVector(conditionalProbability_t aCondProbType);

}
}

/* Special operations:
if(conditionType == Conditions::STEM_TWO_SAMPLES) {
	// Do singleton stuff
}

// unobserved no sampling:
if(conditionType == Conditions::ROOT_SURVIVAL && conditionType == Conditions::ROOT_MRCA &&
   conditionType == Conditions::STEM_SURVIVAL && conditionType == Conditions::STEM_TWO_EXT_SAMPLES) {
   // do unobserved no sampling
}

// singleton no sampling:
if(conditionType == Conditions::STEM_TWO_EXT_SAMPLES) {
	// do singleton no sampling
}*/



#endif /* LIKELIHOOD_CONDITIONTYPES_CONDITIONTYPE_H_ */
