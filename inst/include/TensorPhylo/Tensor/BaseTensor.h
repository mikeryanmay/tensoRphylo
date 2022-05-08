/*
 * BaseTensor.h
 *
 *  Created on: Aug 28, 2019
 *      Author: xaviermeyer
 */

#ifndef TENSOR_BASETENSOR_H_
#define TENSOR_BASETENSOR_H_

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <assert.h>

#include "IncFwdTensor.h"

namespace Tensor {

typedef std::vector<Eigen::MatrixXd> tensor_t;
typedef std::vector< Eigen::SparseMatrix< double, Eigen::RowMajor > >  sparseTensor_t;

typedef std::map< std::vector<unsigned>, double > rbEventMap_t;
typedef rbEventMap_t::const_iterator cItRbEventMap_t;

class BaseTensor {
public:
	BaseTensor();
	virtual ~BaseTensor();

	bool isContantThroughTime() const;
	std::vector<size_t> getDimensions();
	virtual tensor_t& getTensor(double aT) = 0;

	virtual bool isSparse() const = 0;
	virtual sparseTensor_t& getSparseTensor(double aT);

	// reinit the tensor using the referenced elements;
	void update();

	static const double NONZERO_COEFF_SPARSE_THRESHOLD;

protected:


	bool isTimeDependent;
	std::vector<size_t> dimensions;

	void init(); // Init must be called in derived class constructor
	virtual void doResetTensor() = 0;
	virtual void doDefineIfTimeDependent() = 0;
	virtual void doDefineDimensions() = 0;
	virtual void doDefineTensor() = 0;

	virtual double getNonZeroCoefficientRatio() const = 0;

	size_t getEpochIndex(double currentTime, const std::vector<double> &aTimes);

private:

};

} /* namespace Tensor */

#endif /* TENSOR_BASETENSOR_H_ */
