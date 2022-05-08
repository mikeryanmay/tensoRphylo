/*
 * ProbeState.cpp
 *
 *  Created on: Sep 17, 2019
 *      Author: xaviermeyer
 */

#include "ProbeState.h"

#include <sstream>

namespace Likelihood {
namespace Monitor {

ProbeState::ProbeState() {
	hasNegativeImgEIGVal = false;
	time = 0.;
	stiffnessRatio = 0.;
}

ProbeState::~ProbeState() {
}


std::string ProbeState::toString() const {
	std::stringstream ss;

	ss << "Probe at time t= " << time << std::endl;

	ss << "Vector U:" << std::endl;
	for(size_t i=0; i<(size_t)u.size(); ++i) {
		ss << u(i);
		if(i != (size_t)u.size()-1) ss << " ";
	}
	ss << std::endl;

	for(size_t j=0; j<vecP.size(); ++j) {

		ss << "Vector P for edge [" << vecIdEdge[j] << "]" << std::endl;
		for(size_t i=0; i<(size_t)u.size(); ++i) {
			ss << vecP[j](i);
			if(i != (size_t)vecP[j].size()-1) ss << " ";
		}
		ss << " :: scaling = " << scalingFactor[j] << std::endl;
	}

	return ss.str();
}


std::string ProbeState::toDumpFormat() const {
	std::stringstream ss;

	ss << "time" << std::endl;
	ss << std::scientific << time << std::endl;

	ss << "nStates" << std::endl;
	ss << u.size() << std::endl;

	ss << "nEdges" << std::endl;
	ss << vecP.size() << std::endl;

	ss << "U" << std::endl;
	for(size_t i=0; i<(size_t)u.size(); ++i) {
		ss << u(i) << " ";
	}
	ss << std::endl;

	ss << "P" << std::endl;
	for(size_t j=0; j<vecP.size(); ++j) {
		ss << vecIdEdge[j] << " ";
		for(size_t i=0; i<(size_t)u.size(); ++i) {
			ss << vecP[j](i) << " ";
		}
		ss << scalingFactor[j] << std::endl;
	}

	ss << "SR" << std::endl;
	ss << stiffnessRatio << " " << (hasNegativeImgEIGVal ? 1 : 0) << std::endl;

	return ss.str();
}

bool ProbeState::operator<(const ProbeState& other) const {
	return time < other.time;
}

} /* namespace Monitor */
} /* namespace Likelihood */

