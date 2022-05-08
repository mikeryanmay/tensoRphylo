/*
 * Utils.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: xaviermeyer
 */

#include "Utils.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace Likelihood {
namespace StateType {

static const double MIN_SCALING_THRESHOLD = 1.e-3;
static const double MAX_SCALING_THRESHOLD = 1.e+3;
static const double LOG_OF_2 = log(2.);

double rescaleProbabilityVector(Eigen::MatrixXd::ColXpr colProb){
	double scaler = 0.;

	double minCoeff = colProb.minCoeff();
	double maxCoeff = colProb.maxCoeff();

	assert(!(maxCoeff < MIN_SCALING_THRESHOLD && minCoeff > MAX_SCALING_THRESHOLD));

	if(maxCoeff < MIN_SCALING_THRESHOLD) {
		int exponent;
		frexp(maxCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		colProb /= ldexp((double)1., exponent);
	}

	if(minCoeff > MAX_SCALING_THRESHOLD) {
		int exponent;
		frexp(minCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		colProb /= ldexp((double)1., exponent);
	}

	return scaler;
}

double rescaleProbabilityVector(Eigen::Ref<Eigen::MatrixXd>::ColXpr colProb){
	double scaler = 0.;

	double minCoeff = colProb.minCoeff();
	double maxCoeff = colProb.maxCoeff();

	assert(!(maxCoeff < MIN_SCALING_THRESHOLD && minCoeff > MAX_SCALING_THRESHOLD));

	if(maxCoeff < MIN_SCALING_THRESHOLD) {
		int exponent;
		frexp(maxCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		colProb /= ldexp((double)1., exponent);
	}

	if(minCoeff > MAX_SCALING_THRESHOLD) {
		int exponent;
		frexp(minCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		colProb /= ldexp((double)1., exponent);
	}

	return scaler;
}

double rescaleProbabilityVector(Eigen::VectorXd& vecProb) {
	double scaler = 0.;

	double minCoeff = vecProb.minCoeff();
	double maxCoeff = vecProb.maxCoeff();

	assert(!(maxCoeff < MIN_SCALING_THRESHOLD && minCoeff > MAX_SCALING_THRESHOLD));

	if(maxCoeff < MIN_SCALING_THRESHOLD) {
		int exponent;
		frexp(maxCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		vecProb /= ldexp((double)1., exponent);
	}

	if(minCoeff > MAX_SCALING_THRESHOLD) {
		int exponent;
		frexp(minCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		vecProb /= ldexp((double)1., exponent);
	}

	return scaler;
}

double rescaleProbabilityVector(Eigen::Ref<Eigen::VectorXd> vecProb) {
	double scaler = 0.;

	double minCoeff = vecProb.minCoeff();
	double maxCoeff = vecProb.maxCoeff();

	assert(!(maxCoeff < MIN_SCALING_THRESHOLD && minCoeff > MAX_SCALING_THRESHOLD));

	if(maxCoeff < MIN_SCALING_THRESHOLD) {
		int exponent;
		frexp(maxCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		vecProb /= ldexp((double)1., exponent);
	}

	if(minCoeff > MAX_SCALING_THRESHOLD) {
		int exponent;
		frexp(minCoeff, &exponent);
		scaler += exponent*LOG_OF_2;
		vecProb /= ldexp((double)1., exponent);
	}

	return scaler;
}


} /* namespace StateType */
} /* namespace Likelihood */
