#include <Rcpp.h>
#include "DataReader.h"
// #include "Data/Reader/Nexus/NexusParser.h"


///////////////////////
// data reader class //
///////////////////////

using namespace Rcpp;

RCPP_MODULE(DataReaderMod) {

  function("readDelimitedData", &DataReader::readDelimitedData);
  function("readNexusData", &DataReader::readNexusData);

}









//
