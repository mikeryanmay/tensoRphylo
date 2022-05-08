/*
 * NexusParser.h
 *
 *  Created on: Aug 27, 2019
 *      Author: meyerx
 */

#ifndef DATA_READER_NEXUS_NEXUSPARSER_H_
#define DATA_READER_NEXUS_NEXUSPARSER_H_

#include <boost/smart_ptr/shared_ptr.hpp>
#include "ncl/ncl.h"
#include <vector>

namespace Phylogeny {
namespace TSVReader {
class TSVParser;
}

namespace Data {
	class Container;
	typedef boost::shared_ptr<Container> ContainerSharedPtr;
}
}

namespace Phylogeny {
namespace NexusReader {

class NexusParser {
public:
	NexusParser(const std::string &aFileName);
	NexusParser(const std::string &aNexusFileName, const std::string &aTSVFileName);
	~NexusParser();

	bool hasData() const;
	const NxsSimpleTree *getNxsTree() const;
	std::vector<std::string> getTaxaNames() const;
	std::map<int, std::string> getTaxaIdToNameMap() const;
	const std::vector<double>& getProbForTaxaId(int idTaxa) const;
	size_t getNState() const;

	Phylogeny::Data::ContainerSharedPtr createDataContainer() const;

	/* Ugly hack for testing purpose */
	NexusParser(size_t aNState, std::vector<double> &aTimes, std::vector< std::vector<double> > &aInitState);
	bool isSingleBranchHack() const;
	const std::vector<double>& getNodeTimes() const;

private:

	/* Hack */
	bool singleBranchHack, dataFound;
	size_t nState;
	std::vector<double> nodeTimes;

	/* Normal */
	MultiFormatReader nexusReader;
	NxsSimpleTree *nxsTree;

	Phylogeny::TSVReader::TSVParser* tsvParser;

	typedef struct  {
		int idTaxa;
		std::vector<double> initStatesProbs;
	} stateCache_t;

	mutable std::vector<stateCache_t> vecStateCache;

	void init(const std::string &aFileName);
	void init(const std::string &aNexusFileName, const std::string &aTSVFileName);
};

} /* namespace NexusReader */
} /* namespace Phylogeny */

#endif /* DATA_READER_NEXUS_NEXUSPARSER_H_ */
