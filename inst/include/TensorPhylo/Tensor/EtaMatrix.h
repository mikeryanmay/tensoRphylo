/*
 * EtaMatrix.h
 *
 *  Created on: Apr 2, 2020
 *      Author: meyerx
 */

#ifndef TENSOR_ETAMATRIX_H_
#define TENSOR_ETAMATRIX_H_

#include "IncFwdTensor.h"

namespace Tensor {

class EtaMatrix {
public:

	EtaMatrix();
	EtaMatrix(etaStructure_t aEtaStructure);
	virtual ~EtaMatrix();

	void setStructureType(etaStructure_t aEtaStructure);
	etaStructure_t getStructureType() const;

protected:

	etaStructure_t etaStructure;

	bool isQuasseMatrix(const Eigen::MatrixXd &aOriginalMatrix) const;

};

} /* namespace Tensor */

#endif /* TENSOR_ETAMATRIX_H_ */
