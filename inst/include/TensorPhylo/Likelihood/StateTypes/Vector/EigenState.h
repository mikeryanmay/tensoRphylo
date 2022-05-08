/*
 * EigenState.h
 *
 *  Created on: Apr 15, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_STATETYPES_VECTOR_EIGENSTATE_H_
#define LIKELIHOOD_STATETYPES_VECTOR_EIGENSTATE_H_

#include <vector>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/operators.hpp>
#include <boost/numeric/odeint.hpp>

#include "Utils/MemoryPool/EigenCPU.h"
#include "Likelihood/Monitor/ProbeState.h"
#include "Likelihood/ConditionTypes/ConditionType.h"

/*********************** CLASS DEFINITION **********************/

namespace Likelihood {
namespace StateType {
namespace Vector {

class EigenState : public
		boost::additive1< EigenState ,
		boost::additive2< EigenState , double ,
		boost::multiplicative2< EigenState , double > > > {
public:
	EigenState();
	EigenState(const EigenState &aEigenState);
	~EigenState();

	size_t size() const;
	void resize();

	void allocateVecProb();
	void removeVecProb();

	Eigen::VectorXd& getStateProb();
	const Eigen::VectorXd& getStateProb() const;

	void setScaling(double aScaling);
	double getScaling() const;

	void setEdgeMapping(long int aEdgeId);
	long int getEdgeMapping() const;

	// Boost ODEINT algebra
	EigenState& operator+=( const double &val );
	EigenState& operator+=( const EigenState &otherState );
	EigenState& operator*=( const double a );
	EigenState& operator=( const EigenState &otherState );
	double defineNormInf() const;

	// Special operations
	void initMult(double factor, const EigenState &otherState);
	void addMult(double factor, const EigenState &otherState );
	void odeIntRelativeError(double m_eps_abs, double m_eps_rel, double m_a_x, double m_a_dxdt, const EigenState &otherState1, const EigenState &otherState2);

	// Safety step
	void roundNegativeProbabilityToZero();

	void rescaleProbabilities();
	void copyScalingFactors(const EigenState &otherState);

	std::string toString() const;
	Likelihood::Monitor::ProbeState toProbe() const;

private:

	static size_t idSeq;

	const size_t id;

	long int edgeMapping;
	double scaling;

	EigenVectorAlloc* probVec;

	void initMemory();

	friend EigenState operator/( const EigenState &p1 , const EigenState &p2 );
	friend EigenState abs( const EigenState &p );

};

} /* namespace Default */
} /* namespace StateType */
} /* namespace Likelihood */

/*********************** BOOST ODEINT: algebra **********************/

namespace Likelihood {
namespace StateType {
namespace Vector {

EigenState operator/( const EigenState &p1 , const EigenState &p2 );
EigenState abs( const EigenState &p );

} /* namespace Default */
} /* namespace StateType */
} /* namespace Likelihood */

namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< Likelihood::StateType::Vector::EigenState >
{
    typedef double result_type;
    double operator()( const Likelihood::StateType::Vector::EigenState &aState ) const {
    	return aState.defineNormInf();
    }
};

} } }
//]

/*********************** BOOST ODEINT : resize traits **********************/

namespace boost { namespace numeric { namespace odeint {

template< >
struct is_resizeable< Likelihood::StateType::Vector::EigenState > { // declare resizeability
    typedef boost::true_type type;
    const static bool value = type::value;
};

template< >
struct same_size_impl< Likelihood::StateType::Vector::EigenState , Likelihood::StateType::Vector::EigenState > { // define how to check size
    static bool same_size( const Likelihood::StateType::Vector::EigenState &otherState1,
                          				  const Likelihood::StateType::Vector::EigenState &otherState2 )
    {
        return otherState1.size() == otherState2.size();
    }
};

template< >
struct resize_impl< Likelihood::StateType::Vector::EigenState , Likelihood::StateType::Vector::EigenState > { // define how to resize
    static void resize( Likelihood::StateType::Vector::EigenState &otherState1,
                        const Likelihood::StateType::Vector::EigenState &otherState2 )
    {
    	otherState1.resize();
    }
};

} } }


#endif /* LIKELIHOOD_STATETYPES_VECTOR_EIGENSTATE_H_ */
