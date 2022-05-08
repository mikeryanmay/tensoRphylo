/*
 * TimeHeterogenousMatrix.h
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#ifndef TENSOR_TIMEHETEROGENOUSMATRIX_H_
#define TENSOR_TIMEHETEROGENOUSMATRIX_H_

#include <Eigen/Core>
#include "EtaMatrix.h"
#include "BaseTensor.h"

namespace Tensor {

class TimeHeterogenousMatrix: public BaseTensor, public EtaMatrix {
public:
	TimeHeterogenousMatrix(const std::vector<double> &aTimes, const std::vector< Eigen::MatrixXd > &aParameters);
	~TimeHeterogenousMatrix();

	tensor_t& getTensor(double aT);

	bool isSparse() const;
	sparseTensor_t& getSparseTensor(double aT);

private:

	bool sparse;

	const std::vector<double> &times;
	const std::vector< Eigen::MatrixXd > &parameters;

	std::vector< tensor_t > vecTensors;

	std::vector< sparseTensor_t > vecSparseTensor;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_TIMEHETEROGENOUSMATRIX_H_ */
