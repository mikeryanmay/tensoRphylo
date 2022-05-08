/*
 * TimeHeterogenousVector.h
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#ifndef TENSOR_TIMEHETEROGENOUSVECTOR_H_
#define TENSOR_TIMEHETEROGENOUSVECTOR_H_

#include <Eigen/Core>
#include "BaseTensor.h"

namespace Tensor {

class TimeHeterogenousVector: public BaseTensor {
public:
	TimeHeterogenousVector(const std::vector<double> &aTimes, const std::vector< Eigen::VectorXd > &aParameters);
	~TimeHeterogenousVector();

	tensor_t& getTensor(double aT);

	bool isSparse() const;

private:

	bool sparse;

	const std::vector<double> &times;
	const std::vector< Eigen::VectorXd > &parameters;

	std::vector< tensor_t > vecTensors;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_TIMEHETEROGENOUSVECTOR_H_ */
