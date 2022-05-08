/*
 * ParametersSyncParameterContainer.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef PARAMETERS_SYNCPARAMETERSCONTAINER_H_
#define PARAMETERS_SYNCPARAMETERSCONTAINER_H_

#include <Eigen/Core>
#include <vector>

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"

namespace Parameters {

class SyncParameterContainer {
public:
	SyncParameterContainer();
	~SyncParameterContainer();

public: // ugly but convenient for now

	std::vector<double> massSpeciationTimes, massExtinctionTimes, massSamplingTimes, massDestrSamplingTimes;
	std::vector<Eigen::VectorXd> massSpeciationProb, massExtinctionProb, massSamplingProb, massDestrSamplingProb;
	std::vector<Eigen::MatrixXd> massExtinctionStateChangeProb, dummyStateChangeProb;

	friend class Tensor::Factory;

};

typedef boost::shared_ptr<SyncParameterContainer> SyncParameterContainerSharedPtr;

} /* namespace Parameters */

#endif /* PARAMETERS_SYNCPARAMETERSCONTAINER_H_ */
