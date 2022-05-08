/*
 * Manager.h
 *
 *  Created on: Dec 3, 2019
 *      Author: meyerx
 */

#ifndef UTILS_PARALLEL_MANAGER_H_
#define UTILS_PARALLEL_MANAGER_H_

#include <cstddef>
#include <string>
#include <vector>

#include "../Singleton.h"

// Fake define
/*#if not defined(_OPENMP)
	#define _OPENMP
#endif*/

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace Utils {
namespace Parallel {

typedef struct {
	size_t i, j, nI, nJ;
} block_t;

class Manager : public Utils::Singleton<Manager>{

	friend class Utils::Singleton<Manager>;

public:

	void setNThread(size_t aNThread);
	size_t getNThread() const;

	void setMaxNThread(size_t aMaxNThread);
	size_t getMaxNThread() const;

	std::vector<size_t> defineChunks(size_t nThreads, size_t nElements) const;
	std::vector<size_t> defineOptimalChunks(size_t nAvailableThreads, size_t nElements) const;
	std::vector<block_t> defineOptimalBlocks(size_t nAvailableThreads, size_t nElemM, size_t nElemN) const;

	bool isActiveOpenMP() const;
	bool useOpenMP() const;
	bool useBlockParallelOMP(size_t nElements) const;
	bool useBlockParallelOMP(size_t nElements, size_t nAvailableThreads) const;

	void reportNumThread(std::string event);

private:

	static const size_t MIN_ELEMENTS;

	bool enabledOMP;
	size_t nThread, maxNThread;

private:
	Manager();
	~Manager();
};

} /* namespace Parallel */
} /* namespace Utils */

#endif /* UTILS_PARALLEL_MANAGER_H_ */
