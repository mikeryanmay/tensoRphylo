/*
 * ConditionType.cpp
 *
 *  Created on: Sep 11, 2019
 *      Author: xmeyer
 */

#include "ConditionType.h"
#include <assert.h>
#include <iostream>

namespace Likelihood {
namespace Conditions {

conditionalProbability_t intToConditionalProbabilityType(int myInt) {
	switch (myInt) {
		case 0:
			return TIME;
			break;
		case 1:
			return ROOT_SURVIVAL;
			break;
		case 2:
			return ROOT_MRCA;
			break;
		case 3:
			return ROOT_SAMPLING;
			break;
		case 4:
			return STEM_SURVIVAL;
			break;
		case 5:
			return STEM_ONE_SAMPLE;
			break;
		case 6:
			return STEM_TWO_EXT_SAMPLES;
			break;
		case 7:
			return STEM_TWO_SAMPLES;
			break;
		case 8:
			return ROOT_SAMPLING_AND_MRCA;
			break;
		default:
			assert(false && "Unknown conditionalProbabilityType -- unsupported yet.");
			return TIME;
			break;
	}
}

size_t getNVector(conditionalProbability_t aCondProbType) {

	/*ROOT_SURVIVAL=1, 			// Require unobserved - no sampling 			(u_hat)
	ROOT_MRCA=2, 				// Require unobserved - no sampling 			(u_hat)
	ROOT_SAMPLING=3, 			// Require unobserved 							(u)
	STEM_SURVIVAL=4, 			// Require unobserved - no sampling 			(u_hat)
	STEM_ONE_SAMPLE=5, 			// Require unobserved 							(u)
	STEM_TWO_EXT_SAMPLES=6, 	// Require unobserved, singleton - no sampling 	(u_hat, s_hat)
	STEM_TWO_SAMPLES=7 			// Require unobserved, singleton				(u, s)*/

	if(aCondProbType == Likelihood::Conditions::TIME ||
			aCondProbType == Likelihood::Conditions::ROOT_SAMPLING ||
			aCondProbType == Likelihood::Conditions::ROOT_SAMPLING_AND_MRCA ||
			aCondProbType == Likelihood::Conditions::STEM_ONE_SAMPLE) {
		return 1;
	} else if (aCondProbType == Likelihood::Conditions::ROOT_SURVIVAL ||
			aCondProbType == Likelihood::Conditions::ROOT_MRCA ||
			aCondProbType == Likelihood::Conditions::STEM_SURVIVAL ||
			aCondProbType == Likelihood::Conditions::STEM_TWO_SAMPLES) {
		return 2;
	} else if (aCondProbType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
		return 3;
	} else {
		assert(false && "Unknown conditional probability.");
		return 0;
	}
}

}
}
