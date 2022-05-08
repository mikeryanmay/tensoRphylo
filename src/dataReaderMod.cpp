#include <Rcpp.h>
#include "DataReader.h"
// #include "Data/Reader/Nexus/NexusParser.h"


///////////////////////
// data reader class //
///////////////////////

using namespace Rcpp;

RCPP_MODULE(DataReaderMod) {

  Rcpp::function("readDelimitedData", &DataReader::readDelimitedData);
  Rcpp::function("readNexusData", &DataReader::readNexusData);

}









//
