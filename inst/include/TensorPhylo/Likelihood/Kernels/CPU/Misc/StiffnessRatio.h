/*
 * StiffnessRatio.h
 *
 *  Created on: Apr 9, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_KERNELS_CPU_MISC_STIFFNESSRATIO_H_
#define LIKELIHOOD_KERNELS_CPU_MISC_STIFFNESSRATIO_H_

#include "Tensor/IncFwdTensor.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {

class StiffnessRatio {
public:
	StiffnessRatio(Tensor::ContainerSharedPtr aPtrTensorsContainer);
	~StiffnessRatio();

	std::pair<double, bool> estimateStiffnessRatio(double t, const Eigen::VectorXd &u);

private:
	Tensor::ContainerSharedPtr ptrTensors;

};

} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_MISC_STIFFNESSRATIO_H_ */
