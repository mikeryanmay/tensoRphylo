/*
 * Sitmo11RNG.cpp
 *
 *  Created on: 1 oct. 2013
 *      Author: meyerx
 */


#include "Sitmo11RNG.h"

#include "Utils/RNG/Engine/prng_engine.hpp"
#include <boost/random/beta_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/distributions.hpp>

namespace Utils {
namespace RNG {

Sitmo11RNG::Sitmo11RNG() {
	eng = new sitmo::prng_engine();
	setSeed(12345);
}

Sitmo11RNG::Sitmo11RNG(long unsigned int aSeed) {
	eng = new sitmo::prng_engine();
	setSeed(aSeed);
}

Sitmo11RNG::~Sitmo11RNG() {
	delete eng;
}

int Sitmo11RNG::genUniformInt() const {
	boost::random::uniform_int_distribution<int> distribution;
	return distribution(*eng);
}

int Sitmo11RNG::genUniformInt(const int aMin, const int aMax) const {
	if(aMin == aMax) return aMin;
	boost::random::uniform_int_distribution<int> distribution(aMin, aMax);
	return distribution(*eng);
}

double Sitmo11RNG::genUniformDbl() const {
	boost::random::uniform_real_distribution<double> distribution;
	return distribution(*eng);
}

double Sitmo11RNG::genUniformDbl(const double aMin, const double aMax) const {
	if(aMin == aMax) return aMin;
	boost::random::uniform_real_distribution<double> distribution(aMin, aMax);
	return distribution(*eng);
}

double Sitmo11RNG::genGaussian(const double sigma) const {
	boost::random::normal_distribution<double> distribution(0, sigma);
	return distribution(*eng);
}

double Sitmo11RNG::genNormal(const double mu, const double sigma) const {
	boost::random::normal_distribution<double> distribution(mu, sigma);
	return distribution(*eng);
}


double Sitmo11RNG::genLogNormal(const double mu, const double sigma) const {
	boost::random::lognormal_distribution<> distribution(mu, sigma);
	return distribution(*eng);
}

double Sitmo11RNG::genBinomial(const double p, const unsigned int n) const {
	boost::random::binomial_distribution<> distribution(n, p);
	return distribution(*eng);
}

double Sitmo11RNG::genBeta(const double alpha, const double beta) const {
	boost::random::beta_distribution<> distribution(alpha, beta);
	return distribution(*eng);
}

double Sitmo11RNG::genGamma(const double shape, const double scale) const {
	boost::random::gamma_distribution<> distribution(shape, scale);
	return distribution(*eng);
}

double Sitmo11RNG::genExponential(const double lambda) const {
	boost::random::exponential_distribution<> distribution(lambda);
	return distribution(*eng);
}

std::vector<double> Sitmo11RNG::genDirichlet(const std::vector<double> &alphas) const {
	std::vector<double> result(alphas.size(), 0);
	double sum = 0.;
	for(size_t iA=0; iA<alphas.size(); ++iA) {
		result[iA] = genGamma(alphas[iA], 1.);
		sum += result[iA];
	}

	for(size_t iA=0; iA<alphas.size(); ++iA) {
		result[iA] /= sum;
	}
	return result;
}

size_t Sitmo11RNG::getMultinomial(const std::vector<double> &probs) const {
	// NOTE: this assumes that probs is normalized
	double u = genUniformDbl();
	size_t result = 0;
	for(size_t i = 0; i < probs.size(); ++i) {
		u -= probs[i];
		if (u < 0.0) {
			result = i;
			break;
		}
	}
	return result;
}

size_t Sitmo11RNG::getMultinomial(const Eigen::VectorXd &probs) const {
	// NOTE: this assumes that probs is normalized
	double u = genUniformDbl();
	size_t result;
	for(size_t i = 0; i < probs.size(); ++i) {
		u -= probs(i);
		if (u < 0.0) {
			result = i;
			break;
		}
	}
	return result;
}


void Sitmo11RNG::setSeed(const unsigned long int aSeed){
	seed = aSeed;
	eng->set_key(seed, seed/2, seed/4, seed/8);
	eng->set_counter();
}

unsigned long int Sitmo11RNG::getSeed() const {
	return seed;
}

}
}
