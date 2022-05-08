/*
 * EigenState<conditionalProbType>Operations.hpp
 *
 *  Created on: Nov 23, 2019
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATEOPERATIONS_HPP_
#define LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATEOPERATIONS_HPP_

#include <algorithm>

#include <boost/config.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint/util/unit_helper.hpp>

#include "EigenState.h"

namespace Likelihood {
namespace StateType {
namespace Optimized {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType >
struct EigenStateOperations
{
    template < class Fac1 = double >

    struct scale
    {
        const Fac1 m_alpha1;

        scale( Fac1 alpha1 ) : m_alpha1( alpha1 ) { }

        void operator()( EigenState<conditionalProbType> &t1 ) const
        {
            t1 *= m_alpha1;
        }

        typedef void result_type;
    };

    template <   class Fac1 = double >
    struct scale_sum1
    {
        const Fac1 m_alpha1;

        scale_sum1( Fac1 alpha1 ) : m_alpha1( alpha1 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 ) const
        {
			t1.initMult(m_alpha1, t2);
        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 >
    struct scale_sum2
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum2( Fac1 alpha1 , Fac2 alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3) const
        {
          	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 >
    struct scale_sum3
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;

        scale_sum3( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 ) const
        {
          	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;

        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 >
    struct scale_sum4
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;

        scale_sum4( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5) const
        {

          	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;

        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 >
    struct scale_sum5
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;

        scale_sum5( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 , Fac5 alpha5 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6) const
        {
          	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 >
    struct scale_sum6
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;

        scale_sum6( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 , Fac5 alpha5 , Fac6 alpha6 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ){ }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 ,const EigenState<conditionalProbType> &t7) const
        {

        	//t1.resize(t2.size());

        	t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;

        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 >
    struct scale_sum7
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;

        scale_sum7( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 ) const
        {
        	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;

        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 >
    struct scale_sum8
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;

        scale_sum8( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 ) const
        {
        	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb();;

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };

    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 >
    struct scale_sum9
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;

        scale_sum9( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 , const EigenState<conditionalProbType> &t10 ) const
        {
        	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb() +
								m_alpha9*t10.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };

    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 >
    struct scale_sum10
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;

        scale_sum10( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 , Fac10 alpha10 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 , const EigenState<conditionalProbType> &t10 , const EigenState<conditionalProbType> &t11 ) const
        {
        	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb() +
								m_alpha9*t10.getStateProb() +
								m_alpha10*t11.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };


    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 >
    struct scale_sum11
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;

        scale_sum11( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 , const EigenState<conditionalProbType> &t10 , const EigenState<conditionalProbType> &t11 , const EigenState<conditionalProbType> &t12 ) const
        {
        	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb() +
								m_alpha9*t10.getStateProb() +
								m_alpha10*t11.getStateProb() +
								m_alpha11*t12.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };

    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 , class Fac12 = Fac11 >
    struct scale_sum12
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;
        const Fac12 m_alpha12;

        scale_sum12( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 , Fac12 alpha12 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) , m_alpha12( alpha12 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 , const EigenState<conditionalProbType> &t10 , const EigenState<conditionalProbType> &t11 , const EigenState<conditionalProbType> &t12 , const EigenState<conditionalProbType> &t13 ) const
        {
        	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb() +
								m_alpha9*t10.getStateProb() +
								m_alpha10*t11.getStateProb() +
								m_alpha11*t12.getStateProb() +
								m_alpha12*t13.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };

    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 , class Fac12 = Fac11 , class Fac13 = Fac12 >
    struct scale_sum13
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;
        const Fac12 m_alpha12;
        const Fac13 m_alpha13;

        scale_sum13( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 , Fac12 alpha12 , Fac13 alpha13 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) , m_alpha12( alpha12 ) , m_alpha13( alpha13 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 , const EigenState<conditionalProbType> &t10 , const EigenState<conditionalProbType> &t11 , const EigenState<conditionalProbType> &t12 , const EigenState<conditionalProbType> &t13 , const EigenState<conditionalProbType> &t14 ) const
        {
           	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb() +
								m_alpha9*t10.getStateProb() +
								m_alpha10*t11.getStateProb() +
								m_alpha11*t12.getStateProb() +
								m_alpha12*t13.getStateProb() +
								m_alpha13*t14.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };

    template <  class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 , class Fac12 = Fac11 , class Fac13 = Fac12 , class Fac14 = Fac13 >
    struct scale_sum14
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;
        const Fac12 m_alpha12;
        const Fac13 m_alpha13;
        const Fac14 m_alpha14;

        scale_sum14( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 , Fac12 alpha12 , Fac13 alpha13 , Fac14 alpha14 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) , m_alpha12( alpha12 ) , m_alpha13( alpha13 ) , m_alpha14( alpha14 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t4 , const EigenState<conditionalProbType> &t5 , const EigenState<conditionalProbType> &t6 , const EigenState<conditionalProbType> &t7 , const EigenState<conditionalProbType> &t8 , const EigenState<conditionalProbType> &t9 , const EigenState<conditionalProbType> &t10 , const EigenState<conditionalProbType> &t11 , const EigenState<conditionalProbType> &t12 , const EigenState<conditionalProbType> &t13 , const EigenState<conditionalProbType> &t14 , const EigenState<conditionalProbType> &t15 ) const
        {
           	//t1.resize(t2.size());

			t1.getStateProb() = m_alpha1*t2.getStateProb() +
								m_alpha2*t3.getStateProb() +
								m_alpha3*t4.getStateProb() +
								m_alpha4*t5.getStateProb() +
								m_alpha5*t6.getStateProb() +
								m_alpha6*t7.getStateProb() +
								m_alpha7*t8.getStateProb() +
								m_alpha8*t9.getStateProb() +
								m_alpha9*t10.getStateProb() +
								m_alpha10*t11.getStateProb() +
								m_alpha11*t12.getStateProb() +
								m_alpha12*t13.getStateProb() +
								m_alpha13*t14.getStateProb() +
								m_alpha14*t15.getStateProb();

			t1.vecProbToEdgeMapping = t2.vecProbToEdgeMapping;
			t1.vecScaling = t2.vecScaling;
        }

        typedef void result_type;
    };

    template <  class Fac1 = double , class Fac2 = Fac1 >
    struct scale_sum_swap2
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum_swap2( Fac1 alpha1 , Fac2 alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

        void operator()( EigenState<conditionalProbType> &t1 , EigenState<conditionalProbType> &t2 , const EigenState<conditionalProbType> &t3) const
        {
        	assert(false && "Not yet implemented.");
            const EigenState<conditionalProbType> tmp( t1 );
            t1 = m_alpha1 * t2 + m_alpha2 * t3;
            t2 = tmp;
        }

        typedef void result_type;
    };

    /*
     * for usage in for_each2
     *
     * Works with boost::units by eliminating the unit
     */
    template <  class Fac1 = double >
    struct rel_error
    {
        const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

        rel_error( Fac1 eps_abs , Fac1 eps_rel , Fac1 a_x , Fac1 a_dxdt )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt ) { }


        void operator()( EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 ) const
        {
            t3.odeIntRelativeError(m_eps_abs, m_eps_rel, m_a_x, m_a_dxdt, t1, t2);
        }

        typedef void result_type;
    };

    /*
     * for usage in for_each3
     *
     * used in the controller for the rosenbrock4 method
     *
     * Works with boost::units by eliminating the unit
     */
    template <  class Fac1 = double >
    struct default_rel_error
    {
        const Fac1 m_eps_abs , m_eps_rel ;

        default_rel_error( Fac1 eps_abs , Fac1 eps_rel )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) { }


        /*
         * xerr = xerr / ( eps_abs + eps_rel * max( x , x_old ) )
         */
        void operator()( EigenState<conditionalProbType> &t3 , const EigenState<conditionalProbType> &t1 , const EigenState<conditionalProbType> &t2 ) const
        {
        	assert(false && "Not yet implemented.");
            BOOST_USING_STD_MAX();
            using std::abs;
            using boost::numeric::odeint::get_unit_value;
            using boost::numeric::odeint::set_unit_value;
            Fac1 x1 = abs( get_unit_value( t1 ) ) , x2 = abs( get_unit_value( t2 ) );
            set_unit_value( t3 , abs( get_unit_value( t3 ) ) / ( m_eps_abs + m_eps_rel * max BOOST_PREVENT_MACRO_SUBSTITUTION ( x1 , x2 ) ) );
        }

        typedef void result_type;
    };



    /*
     * for usage in reduce
     */

    template <  class Value >
    struct maximum
    {
        template <  class Fac1 , class Fac2 >
        Value operator()( Fac1 t1 , const Fac2 t2 ) const
        {
        	assert(false && "Not yet implemented.");
            using std::abs;
            Value a1 = abs( get_unit_value( t1 ) ) , a2 = abs( get_unit_value( t2 ) );
            return ( a1 < a2 ) ? a2 : a1 ;
        }

        typedef Value result_type;
    };


    template <  class Fac1 = double >
    struct rel_error_max
    {
        const Fac1 m_eps_abs , m_eps_rel;

        rel_error_max( Fac1 eps_abs , Fac1 eps_rel )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel )
        { }

        EigenState<conditionalProbType> operator()( EigenState<conditionalProbType> r , const EigenState<conditionalProbType> &x_old , const EigenState<conditionalProbType> &x , const EigenState<conditionalProbType> &x_err )
        {
        	assert(false && "Not yet implemented.");
            BOOST_USING_STD_MAX();
            using std::abs;
            using boost::numeric::odeint::get_unit_value;
            using boost::numeric::odeint::set_unit_value;
            EigenState<conditionalProbType> tmp = abs( get_unit_value( x_err ) ) / ( m_eps_abs + m_eps_rel * max BOOST_PREVENT_MACRO_SUBSTITUTION ( abs( x_old ) , abs( x ) ) );
            return max BOOST_PREVENT_MACRO_SUBSTITUTION ( r , tmp );
        }
    };


    template <  class Fac1 = double >
    struct rel_error_max2
    {
        const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

        rel_error_max2( Fac1 eps_abs , Fac1 eps_rel , Fac1 a_x , Fac1 a_dxdt )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
        { }

        EigenState<conditionalProbType> operator()( EigenState<conditionalProbType> r , const EigenState<conditionalProbType> &x_old , const EigenState<conditionalProbType> &/*x*/ , const EigenState<conditionalProbType> &dxdt_old , const EigenState<conditionalProbType> &x_err )
        {
        	assert(false && "Not yet implemented.");
            BOOST_USING_STD_MAX();
            using std::abs;
            using boost::numeric::odeint::get_unit_value;
            using boost::numeric::odeint::set_unit_value;
            EigenState<conditionalProbType> tmp = abs( get_unit_value( x_err ) ) /
                    ( m_eps_abs + m_eps_rel * ( m_a_x * abs( get_unit_value( x_old ) ) + m_a_dxdt * abs( get_unit_value( dxdt_old ) ) ) );
            return max BOOST_PREVENT_MACRO_SUBSTITUTION ( r , tmp );
        }
    };


    template <  class Fac1 = double >
    struct rel_error_l2
    {
        const Fac1 m_eps_abs , m_eps_rel;

        rel_error_l2( Fac1 eps_abs , Fac1 eps_rel )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel )
        { }

        template <  class Res , class T1 , class T2 , class T3 >
        Res operator()( Res r , const T1 &x_old , const T2 &x , const T3 &x_err )
        {
        	assert(false && "Not yet implemented.");
            BOOST_USING_STD_MAX();
            using std::abs;
            Res tmp = abs( get_unit_value( x_err ) ) / ( m_eps_abs + m_eps_rel * max BOOST_PREVENT_MACRO_SUBSTITUTION ( abs( x_old ) , abs( x ) ) );
            return r + tmp * tmp;
        }
    };




    template <  class Fac1 = double >
    struct rel_error_l2_2
    {
        const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

        rel_error_l2_2( Fac1 eps_abs , Fac1 eps_rel , Fac1 a_x , Fac1 a_dxdt )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
        { }

        template <  class Res , class T1 , class T2 , class T3 , class T4 >
        Res operator()( Res r , const T1 &x_old , const T2 &/*x*/ , const T3 &dxdt_old , const T4 &x_err )
        {
        	assert(false && "Not yet implemented.");
            using std::abs;
            Res tmp = abs( get_unit_value( x_err ) ) /
                    ( m_eps_abs + m_eps_rel * ( m_a_x * abs( get_unit_value( x_old ) ) + m_a_dxdt * abs( get_unit_value( dxdt_old ) ) ) );
            return r + tmp * tmp;
        }
    };


};


} /* namespace Optimized */
} /* namespace StateType */
} /* namespace Likelihood */


#endif /* LIKELIHOOD_STATETYPES_OPTIMIZED_EIGENSTATEOPERATIONS_HPP_ */
