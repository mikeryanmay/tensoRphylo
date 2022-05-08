/*
 * Utils.h
 *
 *  Created on: Nov 20, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_STATETYPES_UTILS_H_
#define LIKELIHOOD_STATETYPES_UTILS_H_
#include <Eigen/Core>

namespace Likelihood {
namespace StateType {

double rescaleProbabilityVector(Eigen::MatrixXd::ColXpr colProb);
double rescaleProbabilityVector(Eigen::Ref<Eigen::MatrixXd>::ColXpr colProb);
double rescaleProbabilityVector(Eigen::VectorXd& vecProb);
double rescaleProbabilityVector(Eigen::Ref<Eigen::VectorXd> vecProb);

} /* namespace StateType */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_STATETYPES_UTILS_H_ */
