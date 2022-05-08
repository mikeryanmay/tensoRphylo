/*
 * ParametersContainer.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef PARAMETERS_PARAMETERSCONTAINER_H_
#define PARAMETERS_PARAMETERSCONTAINER_H_

#include <Eigen/Core>
#include <vector>

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"

namespace Parameters {

class Container {
public:
	Container();
	~Container();

	void readFromFile(const std::string &filename);
	void writeToFile(const std::string &filename);

public: // ugly but convenient for now
	// For instance

	bool applyTreeLikCorrection;
	int intLikApproximator, intScheme, condType, nThreads;
	double deltaT;
	Eigen::VectorXd rootPrior;

	Eigen::VectorXd lambda, mu, phi, delta;
	Eigen::MatrixXd eta;

	// Most of the tensor operation are of the order: lambda_i * omega_ijk * p_k * p_j
	// omega is assumed to have dimensions :
	// i : omega[x].rows -> first dimension of Eigen::MatrixXd
	// j : omega[x].cols -> second dimension of Eigen::MatrixXd
	// k : omega.size() -> std::vector dimension
	std::vector<Eigen::MatrixXd> omega;

	std::vector<double> massSpeciationTimes, massExtinctionTimes, massSamplingTimes, massDestrSamplingTimes;
	std::vector<Eigen::VectorXd> massSpeciationProb, massExtinctionProb, massSamplingProb, massDestrSamplingProb;
	std::vector<Eigen::MatrixXd> massExtinctionStateChangeProb;

	std::vector<double> synchMonitoring;

	friend class Tensor::Factory;
	friend class SynchronousEvents::Container;

private:

	std::string toString(const std::vector<double> &aVector);
	std::string toString(const std::vector<Eigen::VectorXd> &aVector2);
	std::string toString(const std::vector<Eigen::MatrixXd> &aTensor);
	std::string toString(const Eigen::VectorXd &aVector);
	std::string toString(const Eigen::MatrixXd &aMatrix);

	void fromFile(std::ifstream &myFile, std::vector<double> &aVector);
	void fromFile(std::ifstream &myFile, std::vector<Eigen::VectorXd> &aVector2);
	void fromFile(std::ifstream &myFile, std::vector<Eigen::MatrixXd> &aTensor);
	void fromFile(std::ifstream &myFile, Eigen::VectorXd &aVector);
	void fromFile(std::ifstream &myFile, Eigen::MatrixXd &aMatrix);

};

typedef boost::shared_ptr<Container> ContainerSharedPtr;

} /* namespace Parameters */

#endif /* PARAMETERS_PARAMETERSCONTAINER_H_ */
