/*
 * TSVParser.h
 *
 *  Created on: Nov 20, 2019
 *      Author: xaviermeyer
 */

#ifndef DATA_READER_TSVPARSER_H_
#define DATA_READER_TSVPARSER_H_

#include <string>
#include <map>
#include <vector>

namespace Phylogeny {
namespace TSVReader {

class TSVParser {
public:
	TSVParser(const std::string &aTSVFileName);
	~TSVParser();

	size_t getNState() const;
	const std::vector<std::string>& getTaxaNames() const;

	const std::map< std::string, std::vector<double> >& getDataMap() const;
	bool getValueForKey(std::string key, std::vector<double> &value) const;

private:

	size_t nState;
	std::vector<std::string> taxaNames;
	std::map< std::string, std::vector<double> > data;

	void readFile(const std::string &aTSVFileName);

};

} /* namespace TSVReader */
} /* namespace Phylogeny */

#endif /* DATA_READER_TSVPARSER_H_ */
