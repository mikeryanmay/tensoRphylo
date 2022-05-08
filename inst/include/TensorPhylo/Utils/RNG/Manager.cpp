/*
 * Manager.cpp
 *
 *  Created on: Dec 10, 2019
 *      Author: xaviermeyer
 */

#include "Manager.h"

#include <sys/time.h>

namespace Utils {
namespace RNG {

Manager::Manager() : seed(1) {
	struct timeval time;
	gettimeofday(&time,NULL);
	seed = (time.tv_sec * 1000) + (time.tv_usec / 1000);
	defaultRNG.setSeed(seed);
}

Manager::~Manager() {
}

void Manager::initialize(long int aSeed, size_t aNRNG) {
	seed = aSeed;

	assert(vecRNG.empty() && "You tried to initialize twice the Utils::RNG::Manager.");

	vecRNG.resize(aNRNG);
	for(size_t iR=0; iR<vecRNG.size(); ++iR) {
		vecRNG[iR].setSeed(seed+iR);
	}
}

Sitmo11RNG& Manager::getRNG(size_t aIdRNG) {
	assert(aIdRNG<vecRNG.size());
	return vecRNG[aIdRNG];

}

void Manager::initializeDefaultRNG(long int aSeed) {
	defaultRNG.setSeed(aSeed);
}

Sitmo11RNG& Manager::getDefaultRNG() {
	return defaultRNG;
}

} /* namespace RNG */
} /* namespace Utils */
