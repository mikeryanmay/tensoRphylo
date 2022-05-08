/**
 * @file OutputManager.cpp
 *
 * @date Sep 11, 2019
 * @author meyerx
 * @brief
 */
#include "OutputManager.h"
namespace Utils {
namespace Output {


OutpoutManager::OutpoutManager() : logFileName("tmpLog"), dumpFileName("dump_probes.txt"),verbosityThreshold(MEDIUM_VERB) {
}

OutpoutManager::~OutpoutManager() {
}

void OutpoutManager::setVerbosityThreshold(verboseLevel_t aVerbLvl) {
	verbosityThreshold = aVerbLvl;
}

verboseLevel_t OutpoutManager::getVerbosityThreshold() {
	return verbosityThreshold;
}

void OutpoutManager::setLogFileName(const std::string& aLFName) {
	logFileName = aLFName;
}

const std::string& OutpoutManager::getLogFileName() const {
	return logFileName;
}

void OutpoutManager::setDumpFileName(const std::string& aDFName) {
	dumpFileName = aDFName;
}

const std::string& OutpoutManager::getDumpFileName() const {
	return dumpFileName;
}

bool OutpoutManager::check(verboseLevel_t aVerbLvl) const {
	return verbosityThreshold >= aVerbLvl;
}

} /* namespace Output */
} /* namespace Utils */
