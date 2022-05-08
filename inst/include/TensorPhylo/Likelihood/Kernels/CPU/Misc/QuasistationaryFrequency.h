/*
 * QuasistationaryFrequency.h
 *
 *  Created on: May 6, 2022
 *      Author: mrmay
 */

#ifndef LIKELIHOOD_KERNELS_CPU_MISC_QUASISTATIONARY_H_
#define LIKELIHOOD_KERNELS_CPU_MISC_QUASISTATIONARY_H_

#include "Tensor/IncFwdTensor.h"

namespace Likelihood {
namespace Kernels {
namespace CPU {

Eigen::VectorXd computeQuasiStationaryFrequency(Tensor::ContainerSharedPtr ptrTensors, double t);


} /* namespace CPU */
} /* namespace Kernels */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_KERNELS_CPU_MISC_QUASISTATIONARY_H_ */
