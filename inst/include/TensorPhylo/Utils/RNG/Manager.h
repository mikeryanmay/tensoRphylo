/*
 * Manager.h
 *
 *  Created on: Dec 10, 2019
 *      Author: xaviermeyer
 */

#ifndef UTILS_RNG_MANAGER_H_
#define UTILS_RNG_MANAGER_H_

#include <vector>

#include "../Singleton.h"
#include "Sitmo11RNG.h"

namespace Utils {
namespace RNG {

class Manager : public Utils::Singleton<Manager> {

public:

	void initialize(long int aSeed, size_t aNRNG);
	Sitmo11RNG& getRNG(size_t aIdRNG);

	void initializeDefaultRNG(long int aSeed);
	Sitmo11RNG& getDefaultRNG();

	// Using it!
	// Utils::RNG::Manager::getInstance()->getRNG(1)->DOMYSTUFF()
	// Sitmo11RNG &myRng = getRNG(1);

	// inside your function: Sitmo11RNG &localRNG = Utils::RNG::Manager::getInstance()->getRNG(1);
    // localRNG.genUniformInt();


private:

	long int seed;
	Sitmo11RNG defaultRNG;
	std::vector< Sitmo11RNG > vecRNG;


private:
	Manager();
	~Manager();
	friend class Utils::Singleton<Manager>;

};

} /* namespace RNG */
} /* namespace Utils */

#endif /* UTILS_RNG_MANAGER_H_ */
