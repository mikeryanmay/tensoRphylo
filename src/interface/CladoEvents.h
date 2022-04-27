#ifndef CLADO_H
#define CLADO_H

#include <Rcpp.h>
#include "Interface/DistributionHandler.h"

using namespace Rcpp;
using namespace TensorPhylo::Interface;

// typedef std::map< std::vector<unsigned>, double > eventMap_t;

RCPP_EXPOSED_CLASS(CladoEvents)
class CladoEvents {

  public:

    CladoEvents(size_t dim_);
    CladoEvents(const CladoEvents &other);
    ~CladoEvents(){};

    const size_t&     getDim() const;
    const eventMap_t& getEvents() const;
    const double&     getEvent(std::vector<unsigned> index);
    void              setEvent(std::vector<unsigned> index, double val);


    void          show();

  private:

    void initNoChange();

    // internal data
    size_t     dim;
    eventMap_t data;

};

RCPP_EXPOSED_CLASS(CladoEventsList)
class CladoEventsList {

  public:

    CladoEventsList(size_t i){};
    ~CladoEventsList(){};

    void addCladoEvents(CladoEvents& events) {
      clado_list.push_back(&events);
    }

    CladoEvents getCladoEvents(size_t i) {
      return *(clado_list.at(i - 1));
    }

    size_t size() const {
      return clado_list.size();
    }

    const CladoEvents& at(size_t i) const {
      return *(clado_list.at(i));
    }

    const std::vector<CladoEvents*>& getCladoList() const {
      return clado_list;
    }

    void show() {
      Rcout << "A cladogenetic event array with elements: " << std::endl;
      if ( clado_list.size() > 0 ) {
        // iterate over elements
        for(std::vector<CladoEvents*>::iterator it = clado_list.begin(); it != clado_list.end(); ++it ) {
          (*it)->show();
        }
      }
    }

  private:

    // the vector
    std::vector<CladoEvents*> clado_list;

};


#endif
