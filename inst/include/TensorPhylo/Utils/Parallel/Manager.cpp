/*
 * Manager.cpp
 *
 *  Created on: Dec 3, 2019
 *      Author: meyerx
 */

#include "Manager.h"

#include <cmath>
#include <iostream>

#include<Eigen/Core>

namespace Utils {
namespace Parallel {

const size_t Manager::MIN_ELEMENTS = 4;

Manager::Manager() : nThread(1), maxNThread(1) {

#if defined(_OPENMP)
	enabledOMP = true;
	omp_set_num_threads(nThread);
#else
	enabledOMP = false;
#endif

}

Manager::~Manager() {
}

bool Manager::isActiveOpenMP() const {
	return enabledOMP;
}

bool Manager::useOpenMP() const {
	return enabledOMP && nThread > 1;
}

bool Manager::useBlockParallelOMP(size_t nElements) const {
	size_t optimalNThread = std::min(nElements/MIN_ELEMENTS, maxNThread);
	return useOpenMP() && optimalNThread > 1;
}

bool Manager::useBlockParallelOMP(size_t nElements, size_t nAvailableThreads) const {
	size_t optimalNThread = std::min(nElements/MIN_ELEMENTS, nAvailableThreads);
	return useOpenMP() && optimalNThread > 1;
}

void Manager::setNThread(size_t aNThread) {
	bool hasChanged = aNThread != nThread;
	nThread = aNThread;
#if defined(_OPENMP)
	if(hasChanged) {
		omp_set_num_threads(nThread);
	}
#endif
}

size_t Manager::getNThread() const {
	return nThread;
}

void Manager::setMaxNThread(size_t aMaxNThread) {
	maxNThread = aMaxNThread;
}

size_t Manager::getMaxNThread() const {
	return maxNThread;
}

void Manager::reportNumThread(std::string event) {

#if defined(_OPENMP)
	#pragma omp single
	{
		 std::cout << "Event " << event << ": number of threads in the team - " << omp_get_num_threads() << std::endl;
	}
#endif

}

std::vector<size_t> Manager::defineChunks(size_t nThreads, size_t nElements) const {
	size_t chunkSize = std::floor((double)nElements/(double)nThreads);
	size_t remainder = nElements-nThreads*chunkSize;
	std::vector<size_t> chunks(1,0);
	for(size_t i=0; i<nThreads; ++i) {
		chunks.push_back(chunkSize+chunks.back());
		if(i<remainder) {
			chunks.back() += 1;
		}
	}
	return chunks;
}

std::vector<size_t> Manager::defineOptimalChunks(size_t nAvailableThreads, size_t nElements) const {

	size_t optimalNThread  = std::min(nElements/MIN_ELEMENTS, nAvailableThreads);

	//std::cout << "Available threads : " << nAvailableThreads << " - nOptimalNThread : " << optimalNThread << std::endl;

	size_t chunkSize = std::floor((double)nElements/(double)optimalNThread);
	size_t remainder = nElements-optimalNThread*chunkSize;

	std::vector<size_t> chunks(1,0);

	for(size_t i=0; i<optimalNThread; ++i) {
		chunks.push_back(chunkSize+chunks.back());
		if(i<remainder) {
			chunks.back() += 1;
		}
	}
	return chunks;
}


std::vector<block_t> Manager::defineOptimalBlocks(size_t nAvailableThreads, size_t nElemM, size_t nElemN) const {

	std::vector<block_t> blocks;
	if(nAvailableThreads <= 1) {
		block_t block;
		block.i = 0; // Account for cumulative offset
		block.nI = nElemM;
		block.j = 0; // Account for cumulative offset
		block.nJ = nElemN;
		blocks.push_back(block);
		return blocks;
	}

	//nAvailableThreads *= 4;

	size_t nI = 1;
	size_t nJ = 1;

	size_t sizeI = nElemM;
	size_t sizeJ = nElemN;

	const size_t MIN_SIZE_I = 64;
	const size_t MIN_SIZE_J = 2;
	const size_t MIN_BLOCK_SIZE = MIN_SIZE_I*MIN_SIZE_J;

	int lastDim = 0;

	// Check for chunk size, and then for proc. availability
	bool okMinSizeI = std::floor((float)nElemM/(nI+1.0)) >= MIN_SIZE_I;
	bool okMinSizeJ = std::floor((float)nElemN/(nJ+1.0)) >= MIN_SIZE_J;
	bool okMinBlockI = std::floor(((float)nElemM/(nI+1.0))*sizeJ) >= MIN_BLOCK_SIZE;
	bool okMinBlockJ = std::floor(((float)nElemN/(nJ+1.0))*sizeI) >= MIN_BLOCK_SIZE;
	bool okProcI = ((nI+1)*nJ) <= nAvailableThreads;
	bool okProcJ = (nI*(nJ+1)) <= nAvailableThreads;

	// While we can subdivide one dimension
	while( (okMinSizeI && okMinBlockI && okProcI) || (okMinSizeJ && okMinBlockJ && okProcJ) ) {

		bool isIBiggerThanJ = std::floor((float)sizeI/(float)sizeJ) > 1;
		bool isJBiggerThanI = std::floor((float)sizeJ/(float)sizeI) > 1;

		// We first check by the biggest dimension
		if(isIBiggerThanJ && (okMinSizeI && okMinBlockI && okProcI) ) { // M > 2*N and enough proc
			nI += 1;
			sizeI = std::floor((float)nElemM/(float)nI);
			lastDim = 0;
		} else if(isJBiggerThanI && (okMinSizeJ && okMinBlockJ && okProcJ)) { // N >= M and enough proc
			nJ += 1;
			sizeJ = std::floor((float)nElemN/(float)nJ);
			lastDim = 1;
		} else { // N ~= M or one dimension is stuck
			// We then check what is possible and we alternate the first dimensino each iteration
			if(lastDim == 1) {
				if((okMinSizeI && okMinBlockI && okProcI) ) { //
					nI += 1;
					sizeI = std::floor((float)nElemM/(float)nI);
					lastDim = 0;
				} else if((okMinSizeJ && okMinBlockJ && okProcJ)) {
					nJ += 1;
					sizeJ = std::floor((float)nElemN/(float)nJ);
					lastDim = 1;
				} else {
					break;
				}
			} else {
				if((okMinSizeJ && okMinBlockJ && okProcJ)) {
					nJ += 1;
					sizeJ = std::floor((float)nElemN/(float)nJ);
					lastDim = 1;
				} else if((okMinSizeI && okMinBlockI && okProcI) ) {
					nI += 1;
					sizeI = std::floor((float)nElemM/(float)nI);
					lastDim = 0;
				} else {
					break;
				}
			}
		}
		// update if its possible
		okMinSizeI = std::floor((float)nElemM/(nI+1.0)) >= MIN_SIZE_I;
		okMinSizeJ =  std::floor((float)nElemN/(nJ+1.0)) >= MIN_SIZE_J;
		okMinBlockI = std::floor(((float)nElemM/(nI+1.0))*sizeJ) >= MIN_BLOCK_SIZE;
		okMinBlockJ = std::floor(((float)nElemN/(nJ+1.0))*sizeI) >= MIN_BLOCK_SIZE;
		okProcI = ((nI+1)*nJ) <= nAvailableThreads;
		okProcJ = (nI*(nJ+1)) <= nAvailableThreads;


		/*std::cout << "----------------------------------------------" << std::endl;
		std::cout << nElemM << " x " << nElemN << " ==> " << sizeI << " x " << sizeJ << " =====> nThreads = " << nI << " x " << nJ << " = " << nI*nJ << std::endl;
		std::cout << (okMinSizeI ? "y" : "n") << " -- " << (okMinBlockI ? "y" : "n") << " -- " << (okProcI ? "y" : "n") << std::endl;
		std::cout << (okMinSizeJ ? "y" : "n") << " -- " << (okMinBlockJ ? "y" : "n") << " -- " << (okProcJ ? "y" : "n") << std::endl;
		std::cout << "----------------------------------------------" << std::endl;*/
	}

	/*std::cout << nElemM << " x " << nElemN << " ==> " << sizeI << " x " << sizeJ << " =====> nThreads = " << nI << " x " << nJ << " = " << nI*nJ << std::endl;*/

	// Define blocks
	size_t chunkSizeI = std::floor((double)nElemM/(double)nI);
	size_t remainderI = nElemM-nI*chunkSizeI;

	size_t chunkSizeJ = std::floor((double)nElemN/(double)nJ);
	size_t remainderJ = nElemN-nJ*chunkSizeJ;


	for(size_t i=0; i<nI; ++i) {
		for(size_t j=0; j<nJ; ++j) {
			block_t block;
			block.i = i*sizeI + std::min(i,remainderI); // Account for cumulative offset
			block.nI = sizeI;
			if(i<remainderI) {
				block.nI += 1; // Increase size for the remainder
			}

			block.j = j*sizeJ + std::min(j,remainderJ); // Account for cumulative offset
			block.nJ = sizeJ;
			if(j<remainderJ) {
				block.nJ += 1; // Increase size for the remainder
			}

			//std::cout << "Block " << i << ", " << j << " :: " << block.i << " + " << block.nI << " x " << block.j << " + " << block.nJ << std::endl;
			blocks.push_back(block);
		}
	}

	return blocks;


}

} /* namespace Parallel */
} /* namespace Utils */
