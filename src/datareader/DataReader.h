#ifndef DATAREADER_H
#define DATAREADER_H

#include <Rcpp.h>
using namespace Rcpp;
namespace DataReader {

  NumericMatrix readDelimitedData(std::string file, std::string delim, size_t nstates);
  NumericMatrix readNexusData(std::string file);

}




#endif
