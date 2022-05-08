/*
 * EigenCPU.cpp
 *
 *  Created on: Sep 2, 2019
 *      Author: xaviermeyer
 */

#include "EigenCPU.h"

#include <cassert>
#include <iostream>

// #include "../../Test/Catch2/catch.hpp"
#include "Utils/Parallel/Manager.h"

namespace Utils {
namespace MemoryPool {

EigenCPU::EigenCPU() {
	k=0;
	maxStatesVector = 0;
}

EigenCPU::~EigenCPU() {

	for(size_t i=0; i<vectorPool.size(); ++i) {
		delete vectorPool[i];
	}

	for(size_t i=0; i<matrixPool.size(); ++i) {
		delete matrixPool[i];
	}

	vectorPool.clear();
	matrixPool.clear();

	freeVectorPool.clear();
	freeMatrixPool.clear();

}

void EigenCPU::setNCategories(size_t aK) {

	bool change = false;
	if(aK != k) {
		k = aK;
		change = true;
	}

	if(!vectorPool.empty() && change) {
		for(size_t i=0; i<vectorPool.size(); ++i) {
			assert(vectorPool[i]->free && "Changing the dimension of EigenCPU memory pool while vectors are still allocated.");
			delete vectorPool[i];
		}
		freeVectorPool.clear();
		vectorPool.clear();
	}

	if(!matrixPool.empty()&& change) {
		for(size_t i=0; i<matrixPool.size(); ++i) {
			assert(matrixPool[i]->free && "Changing the dimension of EigenCPU memory pool while vectors are still allocated.");
			delete matrixPool[i];
		}
		freeMatrixPool.clear();
		matrixPool.clear();
	}
}


void EigenCPU::setMaxStatesVector(size_t aMaxStatesVector) {

	bool change = false;
	if(aMaxStatesVector != maxStatesVector) {
		maxStatesVector = aMaxStatesVector;
		change = true;
	}

	if(!matrixPool.empty() && change) {
		for(size_t i=0; i<matrixPool.size(); ++i) {
			assert(matrixPool[i]->free && "Changing the dimension of EigenCPU memory pool while vectors are still allocated.");
			delete matrixPool[i];
		}
		freeMatrixPool.clear();
		matrixPool.clear();
	}
}

size_t EigenCPU::getNMaxStatesVector() const {
	return maxStatesVector;
}

EigenVectorAlloc* EigenCPU::allocateVector() {

	assert(k>0);

	if(!freeVectorPool.empty()) {
		size_t idEva = freeVectorPool.back();
		freeVectorPool.pop_back();
		assert(vectorPool[idEva]->free == true);
		vectorPool[idEva]->free = false;
		//std::cout << "allocate : " << vectorPool[idEva]->id << std::endl;
		return vectorPool[idEva];
	}

	// No free vectors
	EigenVectorAlloc *eva = new EigenVectorAlloc;
	eva->free = false;
	eva->id = vectorPool.size();
	eva->vector.resize(k);
	eva->vector.setZero();
	//std::cout << "Create + allocate : " << eva->id << std::endl;

	vectorPool.push_back(eva);

	return vectorPool.back();

}

EigenMatrixAlloc* EigenCPU::allocateMatrix() {

	assert(k>0);
	assert(maxStatesVector>0);

	if(!freeMatrixPool.empty()) {
		size_t idEma = freeMatrixPool.back();
		freeMatrixPool.pop_back();
		assert(matrixPool[idEma]->free == true);
		matrixPool[idEma]->free = false;

		return matrixPool[idEma];
	}

	// No free matrices
	EigenMatrixAlloc *ema = new EigenMatrixAlloc;
	ema->free = false;
	ema->id = matrixPool.size();
	ema->matrix.resize(k,maxStatesVector);
	ema->matrix.setZero();

	matrixPool.push_back(ema);
	return matrixPool.back();

}

void EigenCPU::freeVector(EigenVectorAlloc *aEVA) {
	assert(aEVA->id < vectorPool.size());
	freeVectorPool.push_back(aEVA->id);
	vectorPool[aEVA->id]->free = true; // Set as free;
	vectorPool[aEVA->id]->vector.setZero();
	//std::cout << "Free : " << aEVA->id << std::endl;
}

void EigenCPU::freeMatrix(EigenMatrixAlloc *aEMA) {

	assert(aEMA->id < matrixPool.size());
	freeMatrixPool.push_back(aEMA->id);
	matrixPool[aEMA->id]->free = true; // Set as free;
	matrixPool[aEMA->id]->matrix.setZero();

}

EigenVectorAlloc* EigenCPU::threadSafeAllocateVector() {

#ifdef _OPENMP
	if(!omp_in_parallel()) {
		return allocateVector();
	} else {
		EigenVectorAlloc* tmpPtr;
		#pragma omp critical (thread_safe_memory_pool_vector)
		{
			tmpPtr = allocateVector();
		}
		return tmpPtr;
	}
#else
	return allocateVector();
#endif

}

EigenMatrixAlloc* EigenCPU::threadSafeAllocateMatrix() {
#ifdef _OPENMP
	if(!omp_in_parallel()) {
		return allocateMatrix();
	} else {
		EigenMatrixAlloc* tmpPtr;
		#pragma omp critical (thread_safe_memory_pool_matrix)
		{
			tmpPtr = allocateMatrix();
		}
		return tmpPtr;
	}
#else
	return allocateMatrix();
#endif
}

void EigenCPU::threadSafeFreeVector(EigenVectorAlloc *aEVA) {
#ifdef _OPENMP
	if(!omp_in_parallel()) {
		freeVector(aEVA);
	} else {
		#pragma omp critical (thread_safe_memory_pool_vector)
		{
			freeVector(aEVA);
		}
	}
#else
	freeVector(aEVA);
#endif
}

void EigenCPU::threadSafeFreeMatrix(EigenMatrixAlloc *aEMA) {
#ifdef _OPENMP
	if(!omp_in_parallel()) {
		freeMatrix(aEMA);
	} else {
		#pragma omp critical (thread_safe_memory_pool_matrix)
		{
			freeMatrix(aEMA);
		}
	}
#else
	freeMatrix(aEMA);
#endif
}


} /* namespace MemoryPool */
} /* namespace Utils */
