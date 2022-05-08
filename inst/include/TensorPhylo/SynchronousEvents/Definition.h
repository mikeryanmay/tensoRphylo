/*
 * Definition.h
 *
 *  Created on: Sep 3, 2019
 *      Author: xaviermeyer
 */

#ifndef SYNCHRONOUSEVENT_DEFINITION_H_
#define SYNCHRONOUSEVENT_DEFINITION_H_

#include <Eigen/Core>
#include <vector>

namespace SynchronousEvents {

class Definition {
public:
	Definition(std::vector<double> &aEventTimes, std::vector<Eigen::VectorXd> &aEventProbabilities, std::vector<Eigen::MatrixXd> &aStateChangeProb);
	~Definition();

	bool hasEvents() const;
	bool hasEvent(double t) const;
	bool areStateChangeInvolved() const;

	std::vector<double>& getEventTimes();
	Eigen::VectorXd& getEventProbability(double t);
	Eigen::MatrixXd& getStateChangeProbability(double t);

public: // Ugly but convenient for now
	std::vector<double> &eventTimes;
	std::vector<Eigen::VectorXd> &eventProbabilities;
	std::vector<Eigen::MatrixXd> &stateChangeProb;
};

} /* namespace SynchronousEvents */

#endif /* SYNCHRONOUSEVENT_DEFINITION_H_ */
