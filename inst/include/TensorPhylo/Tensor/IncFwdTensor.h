/*
 * IncFwdTensor.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef TENSOR_INCFWDTENSOR_H_
#define TENSOR_INCFWDTENSOR_H_
#include <boost/smart_ptr/shared_ptr.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <map>
#include <vector>

namespace Tensor {
	typedef std::vector<Eigen::MatrixXd> tensor_t;
	typedef Eigen::SparseMatrix< double, Eigen::RowMajor > eigenSparseMatrix_t;
	typedef std::vector< eigenSparseMatrix_t >  sparseTensor_t;

	typedef std::map< std::vector<unsigned>, double > rbEventMap_t;
	typedef rbEventMap_t::const_iterator cItRbEventMap_t;

	class EtaMatrix;

	class BaseTensor;
	typedef boost::shared_ptr<BaseTensor> BaseTensorSharedPtr;

	class Container;
	typedef boost::shared_ptr<Container> ContainerSharedPtr;

	class Factory;

	typedef enum {ETA_DENSE = 0,
			          ETA_SPARSE = 1,
				        ETA_QUASSE = 2} etaStructure_t;
}



#endif /* TENSOR_INCFWDTENSOR_H_ */
