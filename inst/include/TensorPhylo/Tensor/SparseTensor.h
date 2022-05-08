/*
 * SparseTensor.h
 *
 *  Created on: March 9, 2020
 *      Author: meyerx
 */

#ifndef TENSOR_SPARSETENSOR_H_
#define TENSOR_SPARSETENSOR_H_

#include "BaseTensor.h"

namespace Tensor {

class SparseTensor: public BaseTensor {
public:
	SparseTensor(size_t aN, const rbEventMap_t &aEventMap);
	~SparseTensor();

	tensor_t& getTensor(double aT);

	bool isSparse() const;
	sparseTensor_t& getSparseTensor(double aT);

private:

	bool sparse;

	size_t N;

	 const rbEventMap_t &rbEventMap;

	sparseTensor_t sparseTensor;

	void doResetTensor();
	void doDefineIfTimeDependent();
	void doDefineDimensions();
	void doDefineTensor();

	double getNonZeroCoefficientRatio() const;

};

} /* namespace Tensor */

#endif /* TENSOR_SPARSETENSOR_H_ */
