/*
 * SimpleTensor.h
 *
 *  Created on: Aug 29, 2019
 *      Author: mrmay
 */

#ifndef TENSOR_SIMPLETENSOR_H_
#define TENSOR_SIMPLETENSOR_H_

#include "BaseTensor.h"

namespace Tensor {

class SimpleTensor: public BaseTensor {
public:
	SimpleTensor(std::vector<Eigen::MatrixXd> &aTensor);
	~SimpleTensor();

	tensor_t& getTensor(double aT);

	bool isSparse() const;
	sparseTensor_t& getSparseTensor(double aT);

private:

	bool sparse;

	tensor_t &originalTensor;

	sparseTensor_t sparseTensor;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_SIMPLETENSOR_H_ */
