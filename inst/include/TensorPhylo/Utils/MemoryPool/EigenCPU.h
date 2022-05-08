/*
 * EigenCPU.h
 *
 *  Created on: Sep 2, 2019
 *      Author: xaviermeyer
 */

#ifndef UTILS_MEMORYPOOL_EIGENCPU_H_
#define UTILS_MEMORYPOOL_EIGENCPU_H_

#include <Eigen/Core>
#include <vector>

#include "IncFwdMemoryPool.h"

namespace Utils {
namespace MemoryPool {

class EigenCPU {
public:

	void setNCategories(size_t aK);
	void setMaxStatesVector(size_t aMaxStatesVector);

	size_t getNMaxStatesVector() const;

	EigenVectorAlloc* allocateVector();
	EigenMatrixAlloc* allocateMatrix();

	void freeVector(EigenVectorAlloc* aEVA);
	void freeMatrix(EigenMatrixAlloc* aEMA);

	EigenVectorAlloc* threadSafeAllocateVector();
	EigenMatrixAlloc* threadSafeAllocateMatrix();

	void threadSafeFreeVector(EigenVectorAlloc* aEVA);
	void threadSafeFreeMatrix(EigenMatrixAlloc* aEMA);

	friend EigenCPU& eigenCPU();

protected:

	size_t k, maxStatesVector;
	std::vector<EigenVectorAlloc* > vectorPool;
	std::vector<EigenMatrixAlloc* > matrixPool;
	std::vector<size_t> freeVectorPool, freeMatrixPool;

private:
	// Defined in the body
	EigenCPU();
	~EigenCPU();
	// Not defined to avoid call
	EigenCPU(const EigenCPU&);
	EigenCPU& operator=(const EigenCPU&);

};

inline EigenCPU& eigenCPU() {
    static EigenCPU instance;
    return instance;
}


} /* namespace MemoryPool */
} /* namespace Utils */

#endif /* UTILS_MEMORYPOOL_EIGENCPU_H_ */
