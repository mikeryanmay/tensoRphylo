/*
 * ParametersSerializer.h
 *
 *  Created on: Mar 15, 2020
 *      Author: meyerx
 */

#ifndef UTILS_PARAMETERSSERIALIZER_H_
#define UTILS_PARAMETERSSERIALIZER_H_
#include <string>
#include <vector>

#include <Eigen/Core>

#include "Tensor/IncFwdTensor.h"

namespace Utils {
namespace Serializer {
namespace Parameters {


std::string toString(const std::vector<double> &aVector);
std::string toString(const std::vector<Eigen::VectorXd> &aVector2);
std::string toString(const std::vector<Eigen::MatrixXd> &aTensor);
std::string toString(const Eigen::VectorXd &aVector);
std::string toString(const Eigen::MatrixXd &aMatrix);
std::string toString(const std::vector< Tensor::rbEventMap_t > &aVecSparseTensor);


void fromFile(std::ifstream &myFile, std::vector<double> &aVector);
void fromFile(std::ifstream &myFile, std::vector<Eigen::VectorXd> &aVector2);
void fromFile(std::ifstream &myFile, std::vector<Eigen::MatrixXd> &aTensor);
void fromFile(std::ifstream &myFile, Eigen::VectorXd &aVector);
void fromFile(std::ifstream &myFile, Eigen::MatrixXd &aMatrix);
void fromFile(std::ifstream &myFile, std::vector< Tensor::rbEventMap_t > &aVecSparseTensor);

} /* namespace Parameters */
} /* namespace Serializer */
} /* namespace Utils */

#endif /* UTILS_PARAMETERSSERIALIZER_H_ */
