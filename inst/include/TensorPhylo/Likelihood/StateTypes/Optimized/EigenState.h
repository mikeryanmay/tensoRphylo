/*
 * EigenState.h
 *
 *  Created on: Nov 18, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATE_H_
#define LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATE_H_

#include <vector>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/operators.hpp>
#include <boost/numeric/odeint.hpp>

#include "Likelihood/Monitor/ProbeState.h"
#include "Utils/MemoryPool/EigenCPU.h"
#include "Likelihood/ConditionTypes/ConditionType.h"

/*********************** CLASS DEFINITION **********************/

namespace Likelihood {
namespace StateType {
namespace Optimized {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
class EigenState : public
		boost::additive1< EigenState< conditionalProbType > ,
		boost::additive2< EigenState< conditionalProbType > , double ,
		boost::multiplicative2< EigenState< conditionalProbType > , double > > > {
public:
	EigenState();
	EigenState(const size_t aNEdges);
	EigenState(const EigenState &aEigenState);
	~EigenState();

	size_t size() const;
	void resize(size_t aNEdges);
	size_t getNVector() const;

	size_t allocateVecProbForEdge(size_t idEdge);
	void removeVecProbForEdge(size_t idEdge);
	void updateVecProbForEdgeMapping(size_t fromIdEdge, size_t toIdEdge);


	Eigen::Ref< Eigen::VectorXd > getUnobservedStateProb();
	Eigen::Ref< Eigen::VectorXd > getUnobservedNoSamplingStateProb();
	Eigen::Ref< Eigen::VectorXd > getSingletonStateProb();
	Eigen::Ref< Eigen::VectorXd > getSingletonNoSamplingStateProb();

	Eigen::Ref< Eigen::VectorXd > getObservedStateProb(size_t iVec);
	Eigen::Ref< Eigen::MatrixXd > getObservedStateProb();
	Eigen::Ref< Eigen::MatrixXd > getStateProb();
	Eigen::Ref< Eigen::MatrixXd > getUnobservedAndObservedStateProb();

	const Eigen::Ref< const Eigen::VectorXd > getUnobservedStateProb() const;
	const Eigen::Ref< const Eigen::VectorXd > getUnobservedNoSamplingStateProb() const;
	const Eigen::Ref< const Eigen::VectorXd > getSingletonStateProb() const;
	const Eigen::Ref< const Eigen::VectorXd > getSingletonNoSamplingStateProb() const;

	const Eigen::Ref< const Eigen::MatrixXd > getObservedStateProb() const;
	const Eigen::Ref< const Eigen::MatrixXd > getStateProb() const;
	const Eigen::Ref< const Eigen::MatrixXd > getUnobservedAndObservedStateProb() const;

	size_t findVecProbIdForEdgeId(size_t iE) const;
	void setVecProbToEdgeMapping(std::vector<int> &aMapping);
	std::vector<int>& getVecProbToEdgeMapping();
	const std::vector<int>& getVecProbToEdgeMapping() const;

	// Boost ODEINT algebra
	EigenState& operator+=( const double &val );
	EigenState& operator+=( const EigenState< conditionalProbType > &state );
	EigenState& operator*=( const double a );
	EigenState& operator=( const EigenState< conditionalProbType > &state );
	double defineNormInf() const;

	// Special operations
	void initMult(double factor, const EigenState< conditionalProbType > &state);
	void addMult(double factor, const EigenState< conditionalProbType > &state );
	void odeIntRelativeError(double m_eps_abs, double m_eps_rel, double m_a_x, double m_a_dxdt, const EigenState< conditionalProbType > &state1, const EigenState< conditionalProbType > &state2);

	// Safety step
	void roundNegativeProbabilityToZero();

	// Scaling
	void rescaleAll();
	double getScalingFactorByVecPos(size_t iV) const;
	double getScalingFactorByEdgeId(size_t iE) const;

	void setScalingFactorByVecPos(size_t iV, double aScalingFactor);
	void setScalingFactorByEdgeId(size_t iE, double aScalingFactor);

	std::string toString() const;
	Likelihood::Monitor::ProbeState toProbe() const;


private:

	static size_t idSeq;
	const size_t id, N_VECTOR;

	std::vector<double> vecScaling;

	EigenMatrixAlloc* stateProb;

	std::vector<int> vecProbToEdgeMapping;

	void initMemory();

	size_t getIdxU() const;
	size_t getIdxUHat() const;
	size_t getIdxS() const;
	size_t getIdxSHat() const;

	template< Likelihood::Conditions::conditionalProbability_t CPT >
	friend EigenState< CPT > operator/( const EigenState< CPT > &p1 , const EigenState< CPT > &p2 );

	template< Likelihood::Conditions::conditionalProbability_t CPT >
	friend EigenState< CPT > abs( const EigenState< CPT > &p );

	template< Likelihood::Conditions::conditionalProbability_t CPT >
	friend struct EigenStateOperations;

};

} /* namespace Optimized */
} /* namespace StateType */
} /* namespace Likelihood */

/*********************** BOOST ODEINT: algebra **********************/

namespace Likelihood {
namespace StateType {
namespace Optimized {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType> operator/( const EigenState<conditionalProbType> &p1 , const EigenState<conditionalProbType> &p2 );

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
EigenState<conditionalProbType> abs( const EigenState<conditionalProbType> &p );

} /* namespace StateType */
} /* namespace Likelihood */
}

namespace boost { namespace numeric { namespace odeint {

template< Likelihood::Conditions::conditionalProbability_t conditionalProbType >
struct vector_space_norm_inf< Likelihood::StateType::Optimized::EigenState< conditionalProbType > >
{
    typedef double result_type;
    double operator()( const Likelihood::StateType::Optimized::EigenState< conditionalProbType > &state ) const {
    	return state.defineNormInf();
    }
};

} } }
//]

/*********************** BOOST ODEINT : resize traits **********************/

namespace boost { namespace numeric { namespace odeint {


template< Likelihood::Conditions::conditionalProbability_t conditionalProbType >
struct is_resizeable< Likelihood::StateType::Optimized::EigenState< conditionalProbType > > { // declare resizeability
    typedef boost::true_type type;
    const static bool value = type::value;
};

template< Likelihood::Conditions::conditionalProbability_t conditionalProbType >
struct same_size_impl< Likelihood::StateType::Optimized::EigenState<conditionalProbType> , Likelihood::StateType::Optimized::EigenState<conditionalProbType> > { // define how to check size
    static bool same_size( const Likelihood::StateType::Optimized::EigenState< conditionalProbType > &state1,
                          				  const Likelihood::StateType::Optimized::EigenState< conditionalProbType > &state2 )
    {
        return state1.size() == state2.size();
    }
};

template< Likelihood::Conditions::conditionalProbability_t conditionalProbType >
struct resize_impl< Likelihood::StateType::Optimized::EigenState<conditionalProbType> , Likelihood::StateType::Optimized::EigenState<conditionalProbType> > { // define how to resize
    static void resize( Likelihood::StateType::Optimized::EigenState< conditionalProbType > &state1,
                        const Likelihood::StateType::Optimized::EigenState< conditionalProbType > &state2 )
    {
    	state1.resize( state2.size() );
    }
};

} } }


#endif /* LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATE_H_ */
