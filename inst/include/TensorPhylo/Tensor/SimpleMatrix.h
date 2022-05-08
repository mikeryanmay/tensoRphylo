/*
 * SimpleMatrix.h
 *
 *  Created on: Aug 29, 2019
 *      Author: mrmay
 */

#ifndef TENSOR_SIMPLEMATRIX_H_
#define TENSOR_SIMPLEMATRIX_H_

#include "EtaMatrix.h"
#include "BaseTensor.h"

namespace Tensor {

class SimpleMatrix: public BaseTensor, public EtaMatrix {
public:
	SimpleMatrix(Eigen::MatrixXd &aMatrix);
	~SimpleMatrix();

	tensor_t& getTensor(double aT);
	tensor_t& getTransposedTensor(double aT);

	bool isSparse() const;
	sparseTensor_t& getSparseTensor(double aT);

private:

	bool sparse;

	tensor_t tensor;

	sparseTensor_t sparseTensor;

	Eigen::MatrixXd &originalMatrix;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_SIMPLEMATRIX_H_ */
