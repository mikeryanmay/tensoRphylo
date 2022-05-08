/*
 * TSVParser.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: xaviermeyer
 */

#include "TSVParser.h"

#include <boost/algorithm/string/case_conv.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Phylogeny {
namespace TSVReader {

TSVParser::TSVParser(const std::string &aTSVFileName) {
	readFile(aTSVFileName);
}

TSVParser::~TSVParser() {

}


size_t TSVParser::getNState() const {
	return nState;
}

const std::vector<std::string>& TSVParser::getTaxaNames() const {
	return taxaNames;
}

const std::map< std::string, std::vector<double> >& TSVParser::getDataMap() const {
	return data;
}

bool TSVParser::getValueForKey(std::string key, std::vector<double> &value) const {

	std::map< std::string, std::vector<double> >::const_iterator it = data.find(key);

	if( it != data.end() ) {
		value = it->second;
		return true;
	} else {
		return false;
	}
}

void TSVParser::readFile(const std::string &aTSVFileName) {

	std::ifstream iFile(aTSVFileName.c_str());
	assert(iFile.good());

	bool nStateFound = false;
	std::string line;

	while(std::getline(iFile, line)) {

		// If empty string, next line
		boost::algorithm::trim(line);
		if(line.empty()) continue;

		// Get tokens
		std::vector<std::string> tokens;
		boost::split(tokens, line, boost::is_any_of("\t "));
		assert(tokens.size() == 2 || (nStateFound && tokens.size() == 1 + nState));

		boost::algorithm::to_lower(tokens[0]);
		if(tokens[0].compare("nstate") == 0 || tokens[0].compare("nstates") == 0) {
			nState = std::atoi(tokens[1].c_str());
			nStateFound = true;
		} else {
			assert(nStateFound);
			for(size_t iC=0; iC < tokens[0].size(); ++iC) {
				if(tokens[0][iC] == '_') tokens[0][iC] = ' ';
			}
			assert(tokens[1].size() > 0);
			if(tokens[1][0] == '?' || tokens[1][0] == '-') {
				taxaNames.push_back(tokens[0]);

				std::vector<double> ones(nState, 1.0);
				std::pair< std::string, std::vector<double> > pairKeyVal = std::make_pair(tokens[0], ones);

				data.insert(pairKeyVal);
			} else {
				if(tokens.size() == 2) {
					int state = std::atoi(tokens[1].c_str());
					assert(state >= 0 && state < (long int)nState);

					std::vector<double> vals(nState, 0.);
					vals[state] = 1.0;

					std::pair< std::string, std::vector<double> > pairKeyVal = std::make_pair(tokens[0], vals);
					taxaNames.push_back(tokens[0]);

					data.insert(pairKeyVal);
				} else if(tokens.size() == (1 + nState)) {
					std::vector<double> vals(nState, 0.);
					for(size_t iV=0; iV<nState; ++iV) {

						long double tmpVal = std::stold(tokens[iV+1]);
						if(tmpVal > std::numeric_limits<double>::max()) {
							tmpVal = std::numeric_limits<double>::max();
						} else if(tmpVal < std::numeric_limits<double>::min()) {
							tmpVal = 0.;
						}

						vals[iV] = tmpVal;

						assert(vals[iV] >= 0);
					}

					std::pair< std::string, std::vector<double> > pairKeyVal = std::make_pair(tokens[0], vals);
					taxaNames.push_back(tokens[0]);

					data.insert(pairKeyVal);
				} else {
					assert(false && "TSV files must contains lines consisting of a) str(taxa_name) int(state), or b) str(taxa_name) vector<double>(state_prob)");
				}
			}
		}
	}

	assert(nStateFound && "Provide a 'nState' variable in your TSV file with the total number of states.");
}

} /* namespace TSVReader */
} /* namespace Phylogeny */
