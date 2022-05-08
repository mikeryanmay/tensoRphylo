/*
 * NexusParser.cpp
 *
 *  Created on: Aug 27, 2019
 *      Author: meyerx
 */

#include "NexusParser.h"

#include <boost/algorithm/string/replace.hpp>
#include <cassert>
#include <cstring>

#include "Data/Reader/Container.h"
#include "../TSV/TSVParser.h"

namespace Phylogeny {
namespace NexusReader {

NexusParser::NexusParser(const std::string &aFileName) :
		nexusReader(-1, NxsReader::WARNINGS_TO_STDERR) {

	tsvParser = NULL;
	dataFound = false;
	singleBranchHack = false;
	nexusReader.SetWarningOutputLevel(NxsReader::ILLEGAL_CONTENT_WARNING);

	init(aFileName);
}


NexusParser::NexusParser(const std::string &aNexusFileName, const std::string &aTSVFileName) :
		nexusReader(-1, NxsReader::WARNINGS_TO_STDERR) {
	dataFound = true;
	singleBranchHack = false;
	nexusReader.SetWarningOutputLevel(NxsReader::ILLEGAL_CONTENT_WARNING);

	init(aNexusFileName, aTSVFileName);
}

NexusParser::NexusParser(size_t aNState, std::vector<double> &aTimes, std::vector< std::vector<double> > &aInitState) :
				nexusReader(-1, NxsReader::WARNINGS_TO_STDERR) {

	tsvParser = NULL;
	singleBranchHack = true;
	dataFound = true;
	nState = aNState;

	assert(aTimes.size()-1 == aInitState.size());

	nexusReader.SetWarningOutputLevel(NxsReader::ILLEGAL_CONTENT_WARNING);

	nodeTimes = aTimes;
	for(size_t iS=0; iS<aInitState.size(); ++iS) {
		stateCache_t stateCache;
		stateCache.idTaxa = iS;
		stateCache.initStatesProbs = aInitState[iS];
		vecStateCache.push_back(stateCache);
	}
	nxsTree = NULL;
}

NexusParser::~NexusParser() {
	nexusReader.DeleteBlocksFromFactories();

	if(tsvParser != NULL) {
		delete tsvParser;
	}

	if(nxsTree != NULL) {
		delete nxsTree;
	}
}

bool NexusParser::hasData() const {
	return dataFound;
}


const NxsSimpleTree *NexusParser::getNxsTree() const {
	assert(!singleBranchHack);
	return nxsTree;
}


Phylogeny::Data::ContainerSharedPtr NexusParser::createDataContainer() const {
	Phylogeny::Data::ContainerSharedPtr ptrData( new Phylogeny::Data::Container());

	std::map<size_t, std::string> mapIdToTaxa;
	std::map<std::string, size_t> mapTaxaToId;

	if(isSingleBranchHack()) { // For the single branch hack
		// Using directly the local NexusParse cache initialized by the dedicated single branch constructor
		for(size_t iS=0; iS<vecStateCache.size(); ++iS) {
			std::stringstream ss;
			ss << vecStateCache[iS].idTaxa;
			ptrData->registerTaxaData(ss.str(), vecStateCache[iS].initStatesProbs);

			mapIdToTaxa[iS] = ss.str();
			mapTaxaToId[ss.str()] = iS;
		}
	} else {
		// Using the information present in the nexus file
		std::vector<std::string> taxaNames(getTaxaNames());

		for(size_t iT=0; iT<taxaNames.size(); ++iT) {
			unsigned int iTaxa = nexusReader.GetTaxaBlock(0)->FindTaxon(NxsString(taxaNames[iT].c_str()));
			mapIdToTaxa[iTaxa] = taxaNames[iT];
			mapTaxaToId[taxaNames[iT]] = iTaxa;

			const std::vector<double>& vecInitProbs = getProbForTaxaId(iTaxa);
			ptrData->registerTaxaData(taxaNames[iT], vecInitProbs);
		}
	}

	ptrData->setTaxaMaps(mapTaxaToId, mapIdToTaxa);

	return ptrData;
}

std::vector<std::string> NexusParser::getTaxaNames() const {
	assert(!singleBranchHack);

	if(tsvParser == NULL) {
		std::vector<std::string> tmpNames(nexusReader.GetTaxaBlock(0)->GetAllLabels());
		for(size_t iN=0; iN<tmpNames.size(); ++iN) {
			boost::replace_all(tmpNames[iN], " ", "_");
		}
		return tmpNames;
	} else {
		return tsvParser->getTaxaNames();
	}
}

std::map<int, std::string> NexusParser::getTaxaIdToNameMap() const {
	assert(!singleBranchHack);

	std::map<int, std::string> mapping;
	std::vector< std::string > taxaNames = getTaxaNames();

	for(size_t iT=0; iT<taxaNames.size(); ++iT) {
		std::string nexusName = taxaNames[iT];
		boost::replace_all(nexusName, "_", " ");
		unsigned int iTaxa = nexusReader.GetTaxaBlock(0)->FindTaxon(NxsString(nexusName.c_str()));
		mapping[iTaxa] = taxaNames[iT];
	}

	return mapping;
}

const std::vector<double>& NexusParser::getProbForTaxaId(int idTaxa) const {
	assert(dataFound);

	for(size_t iF=0; iF<vecStateCache.size(); iF++) {
		if(vecStateCache[iF].idTaxa == idTaxa) {
			return vecStateCache[iF].initStatesProbs;
		}
	}
	assert(!singleBranchHack);

	// we didn't find the data
	stateCache_t stateCache;
	stateCache.idTaxa = idTaxa;

	if(tsvParser == NULL) {
		NxsCharactersBlock *cb = nexusReader.GetCharactersBlock(nexusReader.GetTaxaBlock(0), 0);
		assert(cb->GetNumChar() == 1);

		int iTaxon = idTaxa;
		int iChar = 0;

		bool isGap = cb->IsGapState(iTaxon, iChar);
		bool isMissing = cb->IsMissingState(iTaxon, iChar);

		if(!isMissing && !isGap) { // If it's not a missing or gap state then add states
			int nSymbols = strlen(cb->GetSymbols());
			stateCache.initStatesProbs.assign(nSymbols, 0.);
			int nState = cb->GetNumStates(iTaxon, iChar);
			for(int iState=0; iState<nState; ++iState) {
				int state = cb->GetStateIndex(iTaxon, iChar, iState);
				stateCache.initStatesProbs[state] = 1.;
			}
		} else {
			int nState = strlen(cb->GetSymbols());
			stateCache.initStatesProbs.assign(nState, 1.0);
		}

		vecStateCache.push_back(stateCache);
		return vecStateCache.back().initStatesProbs;
	} else {
		NxsString nxsTaxonName = nexusReader.GetTaxaBlock(0)->GetTaxonLabel(idTaxa);
		std::string taxonName(nxsTaxonName);

		std::vector<double> initValues;
		bool found = tsvParser->getValueForKey(taxonName, initValues);
		assert(found);

		assert(initValues.size() == nState);
		stateCache.initStatesProbs = initValues;

		vecStateCache.push_back(stateCache);
		return vecStateCache.back().initStatesProbs;
	}

	assert(false && "Some node had not initial state for the Single Branch hack");
	return vecStateCache.front().initStatesProbs;
}

size_t NexusParser::getNState() const {
	return nState;
}

bool NexusParser::isSingleBranchHack() const {
	return singleBranchHack;
}

const std::vector<double>& NexusParser::getNodeTimes() const {
	return nodeTimes;
}


void NexusParser::init(const std::string &aFileName) {
	assert(!singleBranchHack);
	try {
		nexusReader.ReadFilepath(aFileName.c_str(), MultiFormatReader::NEXUS_FORMAT);
	} catch(...) {
		nexusReader.DeleteBlocksFromFactories();
		throw;
	}

	NxsTreesBlock *treesBlock = nexusReader.GetTreesBlock(nexusReader.GetTaxaBlock(0), 0);

	const NxsFullTreeDescription &tree = treesBlock->GetFullTreeDescription(0);
	nxsTree = new NxsSimpleTree(tree, 0.0, 0, false);

	if(nexusReader.GetNumCharactersBlocks(nexusReader.GetTaxaBlock(0)) == 0) {
		dataFound = false;
	} else {
		NxsCharactersBlock *cb = nexusReader.GetCharactersBlock(nexusReader.GetTaxaBlock(0), 0);
		assert(cb->GetNumChar() == 1);

		std::string tmp(cb->GetSymbols());
		nState = tmp.size();
		dataFound = true;
	}
}

void NexusParser::init(const std::string &aNexusFileName, const std::string &aTSVFileName) {
	assert(!singleBranchHack);
	try {
		nexusReader.ReadFilepath(aNexusFileName.c_str(), MultiFormatReader::NEXUS_FORMAT);
	} catch(...) {
		nexusReader.DeleteBlocksFromFactories();
		throw;
	}

	NxsTreesBlock *treesBlock = nexusReader.GetTreesBlock(nexusReader.GetTaxaBlock(0), 0);
	const NxsFullTreeDescription &tree = treesBlock->GetFullTreeDescription(0);
	nxsTree = new NxsSimpleTree(tree, 0.0, 0, false);

	tsvParser = new Phylogeny::TSVReader::TSVParser(aTSVFileName);
	nState = tsvParser->getNState();
}

} /* namespace NexusReader */
} /* namespace Phylogeny */
