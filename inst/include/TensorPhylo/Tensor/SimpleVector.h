/*
 * SimpleVector.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef TENSOR_SIMPLEVECTOR_H_
#define TENSOR_SIMPLEVECTOR_H_

#include <Eigen/Core>
#include "BaseTensor.h"

namespace Tensor {

class SimpleVector: public BaseTensor {
public:
	SimpleVector(Eigen::VectorXd &aVector);
	~SimpleVector();

	tensor_t& getTensor(double aT);

	bool isSparse() const;

private:

	bool sparse;
	tensor_t tensor;

	Eigen::VectorXd &originalVector;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_SIMPLEVECTOR_H_ */
