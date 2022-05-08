/*
 * TimeHeterogenousTensor.h
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#ifndef TENSOR_TIMEHETEROGENOUSTENSOR_H_
#define TENSOR_TIMEHETEROGENOUSTENSOR_H_

#include <Eigen/Core>
#include "BaseTensor.h"

namespace Tensor {

class TimeHeterogenousTensor: public BaseTensor {

public:
	TimeHeterogenousTensor(size_t aN, const std::vector<double> &aTimes, const std::vector< rbEventMap_t > &aParameters);
	~TimeHeterogenousTensor();

	tensor_t& getTensor(double aT);

	bool isSparse() const;
	sparseTensor_t& getSparseTensor(double aT);

private:

	bool sparse, empty;

	size_t N;
	const std::vector<double> &times;
	const std::vector< rbEventMap_t > &parameters;

	std::vector< sparseTensor_t > vecSparseTensors;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_TIMEHETEROGENOUSTENSOR_H_ */
