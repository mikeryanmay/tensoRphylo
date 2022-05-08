/**
 * @file OutputManager.h
 *
 * @date Sep 11, 2019
 * @author meyerx
 * @brief
 */
#ifndef VERBOSITY_H_
#define VERBOSITY_H_

#include <stddef.h>
#include <string>

namespace Utils {
namespace Output {

enum verboseLevel_t {
	SILENT_VERB = 0, LOW_VERB = 1, MEDIUM_VERB = 2, HIGH_VERB = 3, DUMP_VERB = 4
};

class OutpoutManager {
public:

	void setVerbosityThreshold(verboseLevel_t aVerbLvl);
	verboseLevel_t getVerbosityThreshold();

	void setLogFileName(const std::string& aLFName);
	const std::string& getLogFileName() const;

	void setDumpFileName(const std::string& aLFName);
	const std::string& getDumpFileName() const;

	bool check(verboseLevel_t aVerbLvl) const;

private:

	std::string logFileName, dumpFileName;
	verboseLevel_t verbosityThreshold;

	// Defined in the body
	OutpoutManager();
	~OutpoutManager();
	// Not defined to avoid call
	OutpoutManager(const OutpoutManager&);
	OutpoutManager& operator=(const OutpoutManager&);

	friend OutpoutManager& outputManager();
};

inline OutpoutManager& outputManager() {
    static OutpoutManager instance;
    return instance;
}

} /* namespace Output */
} /* namespace Utils */
#endif /* VERBOSITY_H_ */

