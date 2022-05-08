#ifndef DATAREADER_H
#define DATAREADER_H

#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string/replace.hpp>
#include "ncl/ncl.h"

using namespace Rcpp;
namespace DataReader {

  inline NumericMatrix readDelimitedData(std::string file, std::string delim, size_t nstates) {

    char del;
    if ( delim == "\t" ) {
      del = '\t';
    } else if ( delim == "," ) {
      del = ',';
    } else {
      stop("Invalid delimited. Use '\t' or ','");
    }

    // read the file
    std::ifstream raw(file);
    std::vector<std::string> taxa;
    std::vector<std::string> data;

    // iterate over lines
    std::string line;
    while ( std::getline(raw, line) ) {

      // get the line
      std::stringstream  lineStream(line);

      // get the taxon and state
      std::string taxon, state;
      std::getline(lineStream, taxon, del);
      std::getline(lineStream, state, del);

      // make sure state isn't empty
      if ( state.empty() ) {
        stop("Error reading table. No data found.");
      }

      // push back
      taxa.push_back(taxon);
      data.push_back(state);

    } // end while over lines

    // get the taxon names
    CharacterVector taxon_names = wrap(taxa);
    size_t ntaxa = taxon_names.size();

    std::vector<std::string> state_labels(nstates);
    for(size_t i = 0; i < nstates; ++i) {
      state_labels.at(i) = std::to_string(i);
    }

    NumericMatrix res(ntaxa, nstates);
    rownames(res) = taxon_names;
    colnames(res) = wrap(state_labels);

    for(size_t i = 0; i < ntaxa; ++i) {

      // get a reference to this data
      const std::string& dat = data.at(i);

      // if missing
      bool isGap     = dat == "-";
      bool isMissing = dat == "?";
      if ( isGap || isMissing ) {

        // didn't find the state label, assume it's missing
        for(size_t j = 0; j < nstates; ++j) {
          res(i,j) = 1.0;
        }

      } else {

        // split the string at slashes
        std::stringstream state_ss(dat);
        std::string substate;
        while ( std::getline(state_ss, substate, '/') ) {

          // turn this state into an int
          size_t state_index = std::stoi( substate );

          // place the data
          res(i, state_index) = 1.0;

        } // end loop over delim

      }

    } // end loop over taxa

    return res;

  }

  inline NumericMatrix readNexusData(std::string file) {

    // create the nexus parser
    MultiFormatReader nexusReader(-1, NxsReader::IGNORE_WARNINGS);

    // try tp read
    try {
      nexusReader.ReadFilepath(file.c_str(), MultiFormatReader::NEXUS_FORMAT);
    } catch(...) {
      nexusReader.DeleteBlocksFromFactories();
      throw;
    }

    // get the character blocks
    NxsCharactersBlock *cb = nexusReader.GetCharactersBlock(nexusReader.GetTaxaBlock(0), 0);

    // get the taxa
    std::vector<std::string> taxa = nexusReader.GetTaxaBlock(0)->GetAllLabels();
    size_t ntaxa = taxa.size();

    // convert the taxa names
    // replace spaces with underscores
    CharacterVector taxa_names(ntaxa);
    for(size_t i = 0; i < ntaxa; ++i) {
      boost::replace_all(taxa[i], " ", "_");
      taxa_names.at(i) = taxa[i];
    }

    // get the character data
    std::string symbols = cb->GetSymbols();
    size_t nstates      = symbols.size();

    // convert to character vector
    CharacterVector labels(nstates);
    std::vector<char> state_chars(symbols.begin(), symbols.end());
    for(size_t i = 0; i < nstates; ++i) {
      labels.at(i) = std::string(1, state_chars[i]);
    }

    // create the matrix
    NumericMatrix data(ntaxa, nstates);
    rownames(data) = taxa_names;
    colnames(data) = labels;

    // loop over taxa
    for(size_t i = 0; i < ntaxa; ++i) {

      // check if the data are missing
      bool isGap     = cb->IsGapState(i, 0);
      bool isMissing = cb->IsMissingState(i, 0);

      if ( isGap || isMissing ) {

        // all missing
        for(size_t s = 0; s < nstates; ++s) {
          data(i,s) = 1.0;
        }

      } else {

        // not missing, but may be ambiguous
        int nState = cb->GetNumStates(i, 0);
        for(size_t s = 0; s < nState; ++s) {
          int state = cb->GetStateIndex(i, 0, s);
          data(i,state) = 1.0;
        }

      } // end assign data for taxon

    } // end loop over taxa

    // make sure to clear out the blocks
    nexusReader.DeleteBlocksFromFactories();

    // return the matrix
    return data;

  }

}

#endif
