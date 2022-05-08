/*
 * ProbeState.h
 *
 *  Created on: Sep 17, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_MONITOR_PROBESTATE_H_
#define LIKELIHOOD_MONITOR_PROBESTATE_H_
#include <Eigen/Core>
#include <string>
#include <vector>

namespace Likelihood {
namespace Monitor {

class ProbeState {
public:
	ProbeState();
	 ~ProbeState();

	 std::string toString() const;
	 std::string toDumpFormat() const;

	 bool operator<(const ProbeState& other) const;

public:

	 bool hasNegativeImgEIGVal;
	 double time, stiffnessRatio;
	 Eigen::VectorXd u;
	 std::vector<size_t> vecIdEdge;
	 std::vector<Eigen::VectorXd> vecP;
	 std::vector<double> scalingFactor;

};

} /* namespace Monitor */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_MONITOR_PROBESTATE_H_ */
