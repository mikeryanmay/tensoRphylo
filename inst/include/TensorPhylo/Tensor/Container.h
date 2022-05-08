/*
 * Container.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef TENSOR_CONTAINER_H_
#define TENSOR_CONTAINER_H_

#include <boost/smart_ptr/shared_ptr.hpp>

#include "IncFwdTensor.h"

namespace Tensor {

class Container {
public:
	Container();
	~Container();

	bool isFullyDefined() const;
	bool areDimensionsMatching() const;
	size_t getNumberOfState() const;

	void setLambda(BaseTensorSharedPtr aPtrLambda);
	void setMu(BaseTensorSharedPtr aPtrMu);
	void setPhi(BaseTensorSharedPtr aPtrPhi);
	void setDelta(BaseTensorSharedPtr aPtrDelta);
	void setOmega(BaseTensorSharedPtr aPtrOmega);
	void setEta(BaseTensorSharedPtr aPtrEta);

	BaseTensorSharedPtr getLambda();
	BaseTensorSharedPtr getMu();
	BaseTensorSharedPtr getPhi();
	BaseTensorSharedPtr getDelta();
	BaseTensorSharedPtr getOmega();
	BaseTensorSharedPtr getEta();

	Eigen::MatrixXd& getEigenVecLambda( double t );
	Eigen::MatrixXd& getEigenVecMu( double t );
	Eigen::MatrixXd& getEigenVecPhi( double t );
	Eigen::MatrixXd& getEigenVecDelta( double t );
	tensor_t& getEigenTensorOmega( double t );
	Eigen::MatrixXd& getEigenMatrixEta( double t );

	etaStructure_t getEtaStructureType() const;

	std::string reportTensorInformations() const;

private:

	BaseTensorSharedPtr lambda, mu, phi, delta; // Parameters vectors
	BaseTensorSharedPtr omega; // Tensor
	BaseTensorSharedPtr eta; // Matrix

};

} /* namespace Tensor */

#endif /* TENSOR_CONTAINER_H_ */
