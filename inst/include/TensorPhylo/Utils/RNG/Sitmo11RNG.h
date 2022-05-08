/*
 * Sitmo11RNG.h
 *
 *  Created on: 1 oct. 2013
 *      Author: meyerx
 */

#ifndef SITMO11RNG_H_
#define SITMO11RNG_H_

#include <Eigen/Core>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "Engine/prng_engine.hpp"

namespace sitmo { class prng_engine; }

namespace Utils {
namespace RNG {


class Sitmo11RNG {
public:
	Sitmo11RNG();
	Sitmo11RNG(long unsigned int aSeed);
	~Sitmo11RNG();

	int genUniformInt() const;
	int genUniformInt(const int aMin, const int aMax) const;

	double genUniformDbl() const;
	double genUniformDbl(const double aMin, const double aMax) const;

	double genGaussian(const double sigma) const;
	double genNormal(const double mu, const double sigma) const;
	double genLogNormal(const double mu, const double sigma) const;

	double genBinomial(const double p, const unsigned int n) const;

	double genBeta(const double alpha, const double beta) const;
	double genGamma(const double shape, const double scale) const;

	double genExponential(const double lambda) const;

	std::vector<double> genDirichlet(const std::vector<double> &alphas) const;

	size_t getMultinomial(const std::vector<double> &probs) const;
	size_t getMultinomial(const Eigen::VectorXd &probs) const;

	void setSeed(const unsigned long int aSeed);
	unsigned long int getSeed() const;

private:
	long unsigned int seed;
	sitmo::prng_engine *eng;
};

}
}

#endif /* SITMO11RNG_H_ */
