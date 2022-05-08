
/*
 * IncFwdMemoryPool.h
 *
 *  Created on: Sep 2, 2019
 *      Author: xaviermeyer
 */

#ifndef UTILS_MEMORYPOOL_INCFWDMEMORYPOOL_H_
#define UTILS_MEMORYPOOL_INCFWDMEMORYPOOL_H_

#include <Eigen/Core>


struct EigenVectorAlloc {
	size_t id;
	bool free;
	Eigen::VectorXd vector;
};

struct EigenMatrixAlloc {
	size_t id;
	bool free;
	Eigen::MatrixXd matrix;
};

class CPU;


#endif /* UTILS_MEMORYPOOL_INCFWDMEMORYPOOL_H_ */
