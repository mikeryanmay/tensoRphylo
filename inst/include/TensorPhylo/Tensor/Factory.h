/*
 * Factoy.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef TENSOR_FACTORY_H_
#define TENSOR_FACTORY_H_

#include <boost/smart_ptr/shared_ptr.hpp>

#include "Tensor/IncFwdTensor.h"
#include "Parameters/IncFwdParameterContainer.h"

namespace Tensor {

class Factory {
public:
	Factory();
	~Factory();

	static Tensor::ContainerSharedPtr createContainerWithTimeHomogenousVectors(Parameters::ContainerSharedPtr ptrParameters);
	static Tensor::ContainerSharedPtr createContainerWithTimeHeterogenousVectors(Parameters::AsyncContainerSharedPtr ptrParameters);

private:
};

} /* namespace Tensor */

#endif /* TENSOR_FACTORY_H_ */
