#ifndef CLADO_H
#define CLADO_H

#include <Rcpp.h>
#include "Interface/DistributionHandler.h"

using namespace Rcpp;
using namespace TensorPhylo::Interface;

// typedef std::map< std::vector<unsigned>, double > eventMap_t;

class CladoEvents {

  public:

    CladoEvents(size_t dim_) : dim(dim_), data() {

      // initialize no-change events
      initNoChange();

    }

    ~CladoEvents(){};

    CladoEvents(const CladoEvents &other)  {
      dim  = other.dim;
      data = other.data;
    }

    const size_t& getDim() const {
      return dim;
    }

    const eventMap_t& getEvents() const {
      return data;
    }

    const double& getEvent(std::vector<unsigned> index) {
      eventMap_t::iterator it = data.find(index);
      if ( it == data.end() ) {
        stop("Event not found.");
      }
      return it->second;
    }

    void setEvent(std::vector<unsigned> index, double val) {

      // check that the dimensions of index are appropriate
      if ( index.size() != 3 ) {
        stop("Event index must be of length 3, i.e., [ancestor, left daughter, right daughter]");
      }

      // check the indices are inside the dimensions
      // (remember that we're indexing starting at zero)
      for(std::vector<unsigned>::iterator it = index.begin(); it != index.end(); ++it) {
        if ( *it < 0 || *it > (dim - 1) ) {
          stop("Cladogenetic event indexes [i,j,k] must must conform to the number of states in the data.");
        }
      }

      // don't allow access to the no-change class
      if ( index[0] == index[1] && index[0] == index[2] ) {
        stop("Please do not attempt to access the no-change (non)event.");
      }

      // check value
      if ( val < 0.0 || val > 1.0 ) {
        stop("Cladogenetic event must have value between 0 and 1 (inclusive).");
      }

      // get the no-change event
      unsigned state_index = index[0];
      std::vector<unsigned> nochange_index(3, state_index);
      eventMap_t::iterator jt = data.find(nochange_index);

      // try to find the value first
      eventMap_t::iterator it = data.find(index);
      double delta;
      if ( it != data.end() ) {

        // compute the change in value
        delta = val - it->second; // this is how much the value is increasing

        // update the value
        it->second = val;

      } else {

        // compute the change in value
        delta = val;

        // otherwise, insert the new value
        data.emplace(index, val);

      }

      // decrease the no-change value accordingly
      jt->second -= delta;

    }

    void show() {

      // the header
      Rcout << "A sparse cladogenetic event array with " << data.size() << " events. <" << this << ">" << std::endl;

      if ( data.size() > 0 ) {
        // iterate over elements
        Rcout << "    (ancestral state -> left daughter state, right daughter state) = value" << std::endl;
        for(eventMap_t::iterator it = data.begin(); it != data.end(); ++it) {
          Rcout << "    (" << it->first[0] + 1 << " -> " << it->first[1] + 1 << ", " << it->first[2] + 1 << ") = " << it->second << std::endl;
        }
      }

    }

  private:

    void initNoChange() {

      // clear the map
      data.clear();

      // add a no-change event for each ancestral state
      for(size_t i = 0; i < dim; ++i) {

        // make the empty index
        std::vector<unsigned> index(3, i);

        // emplace the value
        data.emplace(index, 1.0);

      }

    }

    // internal data
    size_t     dim;
    eventMap_t data;

};

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
