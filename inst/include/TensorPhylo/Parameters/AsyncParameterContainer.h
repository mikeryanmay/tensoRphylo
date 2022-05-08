/*
 * AsyncParameterContainer.h
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#ifndef PARAMETERS_ASYNCPARAMETERSCONTAINER_H_
#define PARAMETERS_ASYNCPARAMETERSCONTAINER_H_

#include <Eigen/Core>
#include <vector>

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"

namespace Parameters {

class AsyncParameterContainer {
public:
	AsyncParameterContainer();
	~AsyncParameterContainer();

	friend class Tensor::Factory;

public:

	std::vector<double> timesLambda, timesMu, timesPhi, timesDelta;
	std::vector< Eigen::VectorXd > vecLambda, vecMu, vecPhi, vecDelta;

	std::vector<double> timesEta, timesOmega;
	std::vector<Eigen::MatrixXd > vecEta;

	size_t nState;
	std::vector< Tensor::rbEventMap_t > vecOmega;

};

typedef boost::shared_ptr<AsyncParameterContainer> AsyncParameterContainerSharedPtr;

} /* namespace Parameters */

#endif /* PARAMETERS_ASYNCPARAMETERSCONTAINER_H_ */
