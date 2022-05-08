/*
 * AutoTuningApproximator.cpp
 *
 *  Created on: Apr 21, 2020
 *      Author: meyerx
 */

#include "AutoTuningApproximator.h"

#include <boost/assign/list_of.hpp>
#include <boost/accumulators/statistics/rolling_count.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>

#include "../../Tensor/BaseTensor.h"
#include "../../Tensor/Container.h"
#include "Data/Structure/Tree.h"
#include "Likelihood/Scheduler/BaseScheduler.h"
#include "IncLikelihoodApproximator.h"

namespace Likelihood {
namespace Approximator {

/* Automatically generated variables: see benchmark/benchmarMusse/scripts/fitModel.py */
const size_t AutoTuningApproximator::thresholdNTips = 128;
const size_t AutoTuningApproximator::thresholdNStates = 128;

//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N), (K**2)*(N**2), (N*K**2.5), (N*K**2.75), (N*K**3)])
const std::vector<double> AutoTuningApproximator::coeffsSmall1A_1T = boost::assign::list_of(0.0)(9.67637306727e-06)(0.0)(0.0)(9.36016535302e-09)(0.0)(4.00302586314e-11)(0.0)(1.84651537033e-10)(1.60721616052e-10);
const std::vector<double> AutoTuningApproximator::coeffsLarge1A_1T = boost::assign::list_of(0.0)(0.0)(0.0)(0.0)(1.07527027432e-08)(9.18935532493e-09)(6.23625851788e-11)(2.28821147782e-09)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall2A_1T = boost::assign::list_of(0.0)(1.01819989816e-05)(7.392453448e-06)(3.68036739503e-07)(0.0)(5.78891672699e-09)(0.0)(0.0)(0.0)(1.38125310137e-10);
const std::vector<double> AutoTuningApproximator::coeffsLarge2A_1T = boost::assign::list_of(0.0)(6.64692140316e-06)(7.49121403821e-05)(0.0)(0.0)(8.74492291598e-09)(0.0)(1.47060592734e-09)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall3A_2T = boost::assign::list_of(0.0)(2.44168784733e-05)(5.0433028083e-06)(7.40701934272e-08)(5.12011837988e-09)(7.09347703716e-09)(0.0)(0.0)(0.0)(7.431141872e-11);
const std::vector<double> AutoTuningApproximator::coeffsLarge3A_2T = boost::assign::list_of(0.0)(3.7835968753e-05)(0.0)(0.0)(4.51429916471e-09)(0.0)(4.46069377441e-11)(6.10211121066e-10)(1.50022497882e-10)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall3A_4T = boost::assign::list_of(0.0)(2.9600114326e-05)(5.33518559219e-06)(3.04298564321e-07)(1.09844886996e-09)(8.18628929654e-09)(0.0)(0.0)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge3A_4T = boost::assign::list_of(0.0)(2.8078915824e-05)(0.00013709336176)(0.0)(1.63266452474e-09)(0.0)(5.09280453491e-11)(0.0)(0.0)(1.05299973918e-11);

const std::vector<double> AutoTuningApproximator::coeffsSmall4A_2T = boost::assign::list_of(0.0)(1.50089787453e-05)(1.82150642583e-06)(4.72248170716e-07)(0.0)(0.0)(1.36021888943e-12)(0.0)(0.0)(9.7758154325e-11);
const std::vector<double> AutoTuningApproximator::coeffsLarge4A_2T = boost::assign::list_of(0.0)(1.50937070304e-05)(8.5871221756e-05)(0.0)(0.0)(4.37818732903e-09)(0.0)(8.49389409906e-10)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall4A_4T = boost::assign::list_of(0.0)(1.62041670321e-05)(3.86502761586e-06)(2.87820999837e-07)(0.0)(0.0)(0.0)(0.0)(0.0)(3.68107103409e-11);
const std::vector<double> AutoTuningApproximator::coeffsLarge4A_4T = boost::assign::list_of(0.0)(1.64611135501e-05)(6.78134089515e-05)(0.0)(0.0)(6.39392652868e-10)(0.0)(5.19163362045e-10)(0.0)(0.0);

// {idApprox, nThreads} => [coeffs]
const std::map< std::pair<size_t, size_t>, std::vector<double> > AutoTuningApproximator::coeffsSmall = boost::assign::map_list_of(std::make_pair(1,1), AutoTuningApproximator::coeffsSmall1A_1T)(std::make_pair(2,1), AutoTuningApproximator::coeffsSmall2A_1T)(std::make_pair(3,2), AutoTuningApproximator::coeffsSmall3A_2T)(std::make_pair(3,4), AutoTuningApproximator::coeffsSmall3A_4T)(std::make_pair(4,2), AutoTuningApproximator::coeffsSmall4A_2T)(std::make_pair(4,4), AutoTuningApproximator::coeffsSmall4A_4T);

const std::map< std::pair<size_t, size_t>, std::vector<double> > AutoTuningApproximator::coeffsLarge = boost::assign::map_list_of(std::make_pair(1,1), AutoTuningApproximator::coeffsLarge1A_1T)(std::make_pair(2,1), AutoTuningApproximator::coeffsLarge2A_1T)(std::make_pair(3,2), AutoTuningApproximator::coeffsLarge3A_2T)(std::make_pair(3,4), AutoTuningApproximator::coeffsLarge3A_4T)(std::make_pair(4,2), AutoTuningApproximator::coeffsLarge4A_2T)(std::make_pair(4,4), AutoTuningApproximator::coeffsLarge4A_4T);

/******************************************************************************************/
/* Automatically generated variables: see benchmark/benchmarQuasse/scripts/fitModel.py */
const size_t AutoTuningApproximator::thresholdNTips_quasse = 128;
const size_t AutoTuningApproximator::thresholdNStates_quasse = 64;

//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N)])
const std::vector<double> AutoTuningApproximator::coeffsSmall1A_1T_quasse = boost::assign::list_of(0.0)(4.7118532589e-06)(2.3072608732e-05)(5.58046917106e-07)(1.53897597114e-08)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge1A_1T_quasse = boost::assign::list_of(0.0)(0.0)(0.0)(0.0)(1.72384811685e-08)(1.53714409378e-08);

const std::vector<double> AutoTuningApproximator::coeffsSmall2A_1T_quasse = boost::assign::list_of(0.0)(0.0)(2.9126932451e-05)(8.31001577621e-07)(0.0)(2.63228375027e-09);
const std::vector<double> AutoTuningApproximator::coeffsLarge2A_1T_quasse = boost::assign::list_of(0.0)(0.0)(0.000222142208785)(0.0)(0.0)(7.1763560553e-09);

const std::vector<double> AutoTuningApproximator::coeffsSmall3A_2T_quasse = boost::assign::list_of(0.0)(2.1416685972e-05)(1.96178836638e-05)(3.81739655478e-07)(3.74791694053e-09)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge3A_2T_quasse = boost::assign::list_of(0.0)(0.0)(2.33216036426e-05)(0.0)(7.40726053398e-09)(4.0916826252e-09);

const std::vector<double> AutoTuningApproximator::coeffsSmall4A_2T_quasse = boost::assign::list_of(0.000648838553731)(2.22651945572e-05)(1.57204364873e-05)(1.05943647425e-06)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge4A_2T_quasse = boost::assign::list_of(0.0)(7.16463843461e-06)(0.000184340485797)(0.0)(0.0)(5.42322404222e-09);

// {idApprox, nThreads} => [coeffs]
const std::map< std::pair<size_t, size_t>, std::vector<double> > AutoTuningApproximator::coeffsSmall_quasse = boost::assign::map_list_of(std::make_pair(1,1), AutoTuningApproximator::coeffsSmall1A_1T_quasse)(std::make_pair(2,1), AutoTuningApproximator::coeffsSmall2A_1T_quasse)(std::make_pair(3,2), AutoTuningApproximator::coeffsSmall3A_2T_quasse)(std::make_pair(4,2), AutoTuningApproximator::coeffsSmall4A_2T_quasse);

const std::map< std::pair<size_t, size_t>, std::vector<double> > AutoTuningApproximator::coeffsLarge_quasse = boost::assign::map_list_of(std::make_pair(1,1), AutoTuningApproximator::coeffsLarge1A_1T_quasse)(std::make_pair(2,1), AutoTuningApproximator::coeffsLarge2A_1T_quasse)(std::make_pair(3,2), AutoTuningApproximator::coeffsLarge3A_2T_quasse)(std::make_pair(4,2), AutoTuningApproximator::coeffsLarge4A_2T_quasse);
/******************************************************************************************/

/******************************************************************************************/
/* Automatically generated variables: see benchmark/benchmarChromosse/scripts/fitModel.py */
const size_t AutoTuningApproximator::thresholdNTips_clado = 256;
const size_t AutoTuningApproximator::thresholdNStates_clado = 256;

//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N), (K**2)*(N**2), (N*K**2.5), ((N**2.5)*(K**2.5))])
const std::vector<double> AutoTuningApproximator::coeffsSmall2A_1T_clado = boost::assign::list_of(0.0)(0.0)(0.000196594059367)(9.52582043053e-07)(0.0)(1.84015873237e-07)(0.0)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge2A_1T_clado = boost::assign::list_of(0.0)(0.0)(0.000196594059367)(9.52582043053e-07)(0.0)(1.84015873237e-07)(0.0)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall1A_1T_clado = boost::assign::list_of(0.0)(0.0)(0.000173425945524)(0.0)(2.51088836156e-08)(1.26041928076e-08)(2.6941525777e-10)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge1A_1T_clado = boost::assign::list_of(0.0)(0.0)(0.000173425945524)(0.0)(2.51088836156e-08)(1.26041928076e-08)(2.6941525777e-10)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall4A_2T_clado = boost::assign::list_of(0.0)(3.83206992495e-06)(8.21797741566e-05)(6.80727943762e-07)(0.0)(9.53765527802e-08)(0.0)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge4A_2T_clado = boost::assign::list_of(0.0)(3.83206992495e-06)(8.21797741566e-05)(6.80727943762e-07)(0.0)(9.53765527802e-08)(0.0)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall3A_2T_clado = boost::assign::list_of(0.0)(9.87520030615e-06)(8.58028584449e-05)(0.0)(1.5314728403e-08)(9.63613579807e-09)(1.23993144822e-10)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge3A_2T_clado = boost::assign::list_of(0.0)(9.87520030615e-06)(8.58028584449e-05)(0.0)(1.5314728403e-08)(9.63613579807e-09)(1.23993144822e-10)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall4A_4T_clado = boost::assign::list_of(0.0)(1.3175047353e-06)(0.000149860934401)(1.3114351739e-07)(0.0)(5.33296960728e-08)(0.0)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge4A_4T_clado = boost::assign::list_of(0.0)(1.3175047353e-06)(0.000149860934401)(1.3114351739e-07)(0.0)(5.33296960728e-08)(0.0)(0.0)(0.0);

const std::vector<double> AutoTuningApproximator::coeffsSmall3A_4T_clado = boost::assign::list_of(0.00289729010237)(5.62096374349e-06)(2.88700989207e-05)(1.09933086964e-07)(1.41986387617e-08)(8.29570662126e-09)(3.75640851048e-11)(0.0)(0.0);
const std::vector<double> AutoTuningApproximator::coeffsLarge3A_4T_clado = boost::assign::list_of(0.00289729010237)(5.62096374349e-06)(2.88700989207e-05)(1.09933086964e-07)(1.41986387617e-08)(8.29570662126e-09)(3.75640851048e-11)(0.0)(0.0);

// {idApprox, nThreads} => [coeffs]
const std::map< std::pair<size_t, size_t>, std::vector<double> > AutoTuningApproximator::coeffsSmall_clado = boost::assign::map_list_of(std::make_pair(2,1), AutoTuningApproximator::coeffsSmall2A_1T_clado)(std::make_pair(1,1), AutoTuningApproximator::coeffsSmall1A_1T_clado)(std::make_pair(4,2), AutoTuningApproximator::coeffsSmall4A_2T_clado)(std::make_pair(3,2), AutoTuningApproximator::coeffsSmall3A_2T_clado)(std::make_pair(4,4), AutoTuningApproximator::coeffsSmall4A_4T_clado)(std::make_pair(3,4), AutoTuningApproximator::coeffsSmall3A_4T_clado);

const std::map< std::pair<size_t, size_t>, std::vector<double> > AutoTuningApproximator::coeffsLarge_clado = boost::assign::map_list_of(std::make_pair(2,1), AutoTuningApproximator::coeffsLarge2A_1T_clado)(std::make_pair(1,1), AutoTuningApproximator::coeffsLarge1A_1T_clado)(std::make_pair(4,2), AutoTuningApproximator::coeffsLarge4A_2T_clado)(std::make_pair(3,2), AutoTuningApproximator::coeffsLarge3A_2T_clado)(std::make_pair(4,4), AutoTuningApproximator::coeffsLarge4A_4T_clado)(std::make_pair(3,4), AutoTuningApproximator::coeffsLarge3A_4T_clado);
/******************************************************************************************/


const size_t AutoTuningApproximator::ROLLING_MEAN_WINDOW_SIZE = 100;
const size_t AutoTuningApproximator::AUTO_TUNING_LENGTH = 25;

AutoTuningApproximator::AutoTuningApproximator(
		Likelihood::Integrator::integrationScheme_t aIntScheme,
		Conditions::conditionalProbability_t aConditionType,
		Phylogeny::Data::ContainerSharedPtr aPtrData,
		Scheduler::SchedulerSharedPtr aPtrScheduler,
		SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
		Tensor::ContainerSharedPtr aPtrTensorCont) :
		BaseApproximator(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont) {

	likCounter = 0;
	avgTimeAtSelection = -1;
	currentApproximatorType = predictBestIntegrator();
	ptrApprox = createApproximator(currentApproximatorType);
}

AutoTuningApproximator::~AutoTuningApproximator() {
}

void AutoTuningApproximator::enableTreeLikelihoodCorrection() {
	ptrApprox->enableTreeLikelihoodCorrection();
}

void AutoTuningApproximator::disableTreeLikelihoodCorrection() {
	ptrApprox->disableTreeLikelihoodCorrection();
}

void AutoTuningApproximator::setQuasistationaryFrequencyMode(bool setActive) {
	ptrApprox->setQuasistationaryFrequencyMode(setActive);
}

Eigen::VectorXd AutoTuningApproximator::getRootFrequency(double t) {
	return ptrApprox->getRootFrequency(t);
}

void AutoTuningApproximator::setConditionalProbabilityCompatibilityMode(bool setActive) {
	ptrApprox->setConditionalProbabilityCompatibilityMode(setActive);
}

void AutoTuningApproximator::setPriorStateProbability(const Eigen::VectorXd &aPriorStateProbability) {
	ptrApprox->setPriorStateProbability(aPriorStateProbability);
}

void AutoTuningApproximator::orderProbes() {
	ptrApprox->orderProbes();
}

const std::vector<Likelihood::Monitor::ProbeState>& AutoTuningApproximator::getObservedProbesState() const {
	return ptrApprox->getObservedProbesState();
}

const std::vector<double>& AutoTuningApproximator::getIntegrationTimes() const {
	return ptrApprox->getIntegrationTimes();
}

double AutoTuningApproximator::approximateLogLikelihood() {

	size_t currentNThreads = Utils::Parallel::Manager::getInstance()->getNThread();
	size_t maxNThreads = Utils::Parallel::Manager::getInstance()->getMaxNThread();
	size_t nextNThreadToTry = defineNextSettingToTry(currentNThreads, maxNThreads);

	//std::cout << "NThread = " << currentNThreads << std::endl;
	cp.startTime(0);
	double llik = ptrApprox->approximateLogLikelihood();
	cp.endTime(0);

	if(!nThreadPerf.empty()) {
		if(avgTimeAtSelection < 0.) { // We are tuning - frequent update
			nThreadPerf[nextNThreadToTry](cp.duration(0));
		} else if(likCounter % 20 == 0) { // We are no longer tuning... checking every 20th iterations
			nThreadPerf[nextNThreadToTry](cp.duration(0));
			double ratioTimes = boost::accumulators::rolling_mean(nThreadPerf[currentNThreads]) / avgTimeAtSelection;
			//std::cout << likCounter << " -- " << ratioTimes << std::endl;
			if( avgTimeAtSelection > 0. && (ratioTimes < 1./2. || ratioTimes > 2./1.) ) {
				//std::cout << "[AUTOTUNER][DEBUG] retuning." << std::endl;
				//std::cout << "Reset? : " << avgTimeAtSelection << " vs " << boost::accumulators::rolling_mean(nThreadPerf[currentNThreads]) << std::endl;
				likCounter = 0;
				avgTimeAtSelection = -1.;
				nThreadPerf.assign(Utils::Parallel::Manager::getInstance()->getMaxNThread()+1, accumulator_t(boost::accumulators::tag::rolling_window::window_size = ROLLING_MEAN_WINDOW_SIZE));
			}
		}


	}

	likCounter++;
	//std::cout << "Time = " << cp.duration(0) << std::endl;
	//std::cout << "Lik counter : " << likCounter << std::endl;

	// select the best approximator amon the one tested
#if defined(_OPENMP)
	if((currentApproximatorType == PARALLEL_OPTIMIZED || currentApproximatorType == PARALLEL_BRANCHWISE) && likCounter <= AUTO_TUNING_LENGTH*ROLLING_MEAN_WINDOW_SIZE) {
#else
	if(false) {
#endif
		//std::cout << "Time = " << cp.duration(0) << std::endl;
	    if(boost::accumulators::rolling_count(nThreadPerf[currentNThreads]) >= ROLLING_MEAN_WINDOW_SIZE/2) {
			size_t bestSetting = 1;
			double avgTime = 1e100; //std::numeric_limits<double>::max();
			if(boost::accumulators::rolling_count(nThreadPerf[bestSetting]) >= ROLLING_MEAN_WINDOW_SIZE/2) {
				avgTime = boost::accumulators::rolling_mean(nThreadPerf[bestSetting]);
			}

			for(size_t iT=2; iT<=maxNThreads; ++iT) {
				if(boost::accumulators::rolling_count(nThreadPerf[iT]) < ROLLING_MEAN_WINDOW_SIZE/2) continue;
				// We only change if the gain of performance > 3%
				double expectedSpeedup = (double)iT/(double)bestSetting; // from precedent
				//std::cout << avgTime << " -- " << boost::accumulators::rolling_mean(nThreadPerf[iT]) << std::endl;
				double observedSpeedup = avgTime / boost::accumulators::rolling_mean(nThreadPerf[iT]);
				double percentGainVsExpected = 100.*(observedSpeedup-1.0)/(expectedSpeedup-1.0);
				//std::cout << "Candidate = " << iT << " => " << boost::accumulators::rolling_mean(nThreadPerf[iT]) << " :: expectedSpeedup = " << expectedSpeedup << " - observed = " << observedSpeedup << " - percent = " << percentGainVsExpected << std::endl;
				if(percentGainVsExpected > 10.) { // We see some speedup that is at least 10% of the expected
					bestSetting = iT;
					avgTime = boost::accumulators::rolling_mean(nThreadPerf[iT]);
				}
			}
			//std::cout << "Selecting = " << bestSetting << " => " << boost::accumulators::rolling_mean(nThreadPerf[bestSetting]) << std::endl;
			if(avgTime < 1e100) {
				currentNThreads = bestSetting;
			} // Otherwise leave untouched
	    }
	    Utils::Parallel::Manager::getInstance()->setNThread(currentNThreads);
	    //std::cout << "-----------------------------------------------------------"  << std::endl;
	    //std::cout.flush();

	    if( likCounter == AUTO_TUNING_LENGTH*ROLLING_MEAN_WINDOW_SIZE) {
	    	avgTimeAtSelection = boost::accumulators::rolling_mean(nThreadPerf[currentNThreads]);
	    }
	}

	return llik;
}

void AutoTuningApproximator::setDefaultDeltaT(double aDeltaT) {
	ptrApprox->setDefaultDeltaT(aDeltaT);
}

size_t AutoTuningApproximator::getTotalNumberOfIntegrationSteps() const {
	return ptrApprox->getTotalNumberOfIntegrationSteps();
}


bool AutoTuningApproximator::isLikelihoodDroppingWithNT(size_t currentNThreads, size_t maxNThreads){

	double avgCurTime = boost::accumulators::rolling_mean(nThreadPerf[currentNThreads]);

	// Check if dropping
	size_t cntVisited = 0;
	bool areDropping = true;
	for(size_t iT=currentNThreads+1; iT<=maxNThreads; ++iT) {
		if(boost::accumulators::rolling_count(nThreadPerf[iT]) > ROLLING_MEAN_WINDOW_SIZE/2) {
			cntVisited++;
			areDropping = areDropping && (0.98*boost::accumulators::rolling_mean(nThreadPerf[iT]) > avgCurTime);
		}
	}
	if(cntVisited>0 && areDropping) return areDropping;

	// Check that we don't only test the maxNThread setting
	bool isOnlyMaxThreadTested = currentNThreads == maxNThreads;
	if(!isOnlyMaxThreadTested) {
		return false;
	} else {
		for(size_t iT=1; iT<currentNThreads; ++iT) {
			isOnlyMaxThreadTested = isOnlyMaxThreadTested && boost::accumulators::rolling_count(nThreadPerf[iT]) < ROLLING_MEAN_WINDOW_SIZE/2;
		}
	}
	return isOnlyMaxThreadTested;
}

size_t AutoTuningApproximator::defineNextSettingToTry(size_t currentNThreads, size_t maxNThreads) {
	size_t toTryNThread = currentNThreads;
#if defined(_OPENMP)
	if((currentApproximatorType == PARALLEL_OPTIMIZED || currentApproximatorType == PARALLEL_BRANCHWISE) && likCounter < AUTO_TUNING_LENGTH*ROLLING_MEAN_WINDOW_SIZE) {
#else
	if(false) {
#endif
		/*std::cout << maxNThreads << " -- " << nThreadPerf.size() << std::endl;
		std::cout << "Current NThread = " << currentNThreads << std::endl;
		for(size_t iT=1; iT<=maxNThreads; ++iT) {
			std::cout << iT << " -> " << boost::accumulators::rolling_mean(nThreadPerf[iT]) << "\t";
		}
		std::cout << std::endl;*/

		if(boost::accumulators::rolling_count(nThreadPerf[currentNThreads]) >= ROLLING_MEAN_WINDOW_SIZE/2) {

			bool isDropping = isLikelihoodDroppingWithNT(currentNThreads, maxNThreads);

			if(!isDropping) {
				//std::cout << "Not dropping" << std::endl;
				// Try some other setting action
				size_t action = likCounter % 6;
				if(currentNThreads+1 <= maxNThreads && (action % 3) == 0) { 		// Action 0 and 3
					toTryNThread = currentNThreads+1;
					Utils::Parallel::Manager::getInstance()->setNThread(toTryNThread);
				} else if (currentNThreads+2 <= maxNThreads && (action % 3) == 1) { // Action 1 and 4
					toTryNThread = currentNThreads+2;
					Utils::Parallel::Manager::getInstance()->setNThread(toTryNThread);
				} else { 																									// Action 2 and 5
					// Do nothing - refine the time approx for current setting
				}
			} else {
				//std::cout << "Dropping" << std::endl;
				// Try some other setting action
				size_t action = likCounter % 8;
				if(currentNThreads+1 <= maxNThreads && (action % 4) == 0) { 		// Action 0 and 4
					toTryNThread = currentNThreads+1;
					Utils::Parallel::Manager::getInstance()->setNThread(toTryNThread);
				} else if (currentNThreads+2 <= maxNThreads && (action % 4) == 1 && action == 3) { // Action 1, 3 and 5
					toTryNThread = currentNThreads+2;
					Utils::Parallel::Manager::getInstance()->setNThread(toTryNThread);
				} else if (currentNThreads > 1 && action == 2) { // Action 2
					toTryNThread = currentNThreads-1;
					Utils::Parallel::Manager::getInstance()->setNThread(toTryNThread);
				} else if (currentNThreads > 2 && action == 6) { // Action 6
					toTryNThread = currentNThreads-2;
					Utils::Parallel::Manager::getInstance()->setNThread(toTryNThread);
				} else { 																									// Action 3
					// Do nothing - refine the time approx for current setting
				}
			}
		}
		//std::cout << "Trying : " << toTryNThread << std::endl;
	}
	return toTryNThread;
}

double AutoTuningApproximator::predictLikelihoodEvaluationTimeMusse(size_t nTips, size_t kStates, const std::vector<double> &coeffs) const {
	/* Based on models fitted with least square regression: see benchmark/benchmarMusse/scripts/fitModel.py */
	//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N), (K**2)*(N**2), (N*K**2.5), (N*K**2.75), (N*K**3)])
	assert(coeffs.size() == 10);

	double predictedTime = coeffs[0];
	predictedTime += coeffs[1]*(nTips);
	predictedTime += coeffs[2]*(kStates);
	predictedTime += coeffs[3]*(kStates*nTips);
	predictedTime += coeffs[4]*(kStates*std::pow(nTips, 2.));
	predictedTime += coeffs[5]*(std::pow(kStates, 2.)*nTips);
	predictedTime += coeffs[6]*(std::pow(kStates, 2.)*std::pow(nTips, 2.));
	predictedTime += coeffs[7]*(std::pow(kStates, 2.5)*nTips);
	predictedTime += coeffs[8]*(std::pow(kStates, 2.75)*nTips);
	predictedTime += coeffs[9]*(std::pow(kStates, 3.)*nTips);

	return predictedTime;
}


double AutoTuningApproximator::predictLikelihoodEvaluationTimeClasse(size_t nTips, size_t kStates, const std::vector<double> &coeffs) const {
	/* Based on models fitted with least square regression: see benchmark/benchmarChromosse/scripts/fitModel.py */
	//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N), (K**2)*(N**2), (N*K**2.5), ((N**2.5)*(K**2.5))])
	assert(coeffs.size() == 9);

	double predictedTime = coeffs[0];
	predictedTime += coeffs[1]*(nTips);
	predictedTime += coeffs[2]*(kStates);
	predictedTime += coeffs[3]*(kStates*nTips);
	predictedTime += coeffs[4]*(kStates*std::pow(nTips, 2.));
	predictedTime += coeffs[5]*(std::pow(kStates, 2.)*nTips);
	predictedTime += coeffs[6]*(std::pow(kStates, 2.)*std::pow(nTips, 2.));
	predictedTime += coeffs[7]*(std::pow(kStates, 2.5)*nTips);
	predictedTime += coeffs[8]*(std::pow(kStates, 2.5)*std::pow(nTips, 2.5));

	return predictedTime;
}

double AutoTuningApproximator::predictLikelihoodEvaluationTimeQuasse(size_t nTips, size_t kStates, const std::vector<double> &coeffs) const {
	/* Based on models fitted with least square regression: see benchmark/benchmarChromosse/scripts/fitModel.py */
	//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N)])
	assert(coeffs.size() == 6);

	double predictedTime = coeffs[0];
	predictedTime += coeffs[1]*(nTips);
	predictedTime += coeffs[2]*(kStates);
	predictedTime += coeffs[3]*(kStates*nTips);
	predictedTime += coeffs[4]*(kStates*std::pow(nTips, 2.));
	predictedTime += coeffs[5]*(std::pow(kStates, 2.)*nTips);

	return predictedTime;
}


//coefficients : ([K*0+1, N, K, N*K, (K)*(N**2), (K**2)*(N)])

approximatorVersion_t AutoTuningApproximator::predictBestIntegrator() {

	size_t nTips = ptrScheduler->getPtrTree()->getTerminalNodes().size();
	size_t kStates = ptrTensorCont->getNumberOfState();
	bool isParallel = Utils::Parallel::Manager::getInstance()->isActiveOpenMP() && Utils::Parallel::Manager::getInstance()->getMaxNThread() > 1;

	std::pair<size_t, size_t> bestSetting;
	double minTime = std::numeric_limits<double>::max();

	// dudt: step 2
	if(ptrTensorCont->getEtaStructureType() == Tensor::ETA_QUASSE) {
		/*if(!isParallel || kStates < 64 ) { // Sequential
			bestSetting = std::make_pair(SEQUENTIAL_BRANCHWISE, 1);
		} else { // Parallel
			if(nTips <= 128) {
				bestSetting = std::make_pair(PARALLEL_OPTIMIZED, 4);
			} else {
				bestSetting = std::make_pair(PARALLEL_BRANCHWISE, 4);
			}
		}*/
		typedef std::map< std::pair<size_t, size_t>, std::vector<double> >::const_iterator cItMap;
		if(nTips <= thresholdNTips_quasse && kStates <= thresholdNStates_quasse) { // small-scale model fitted
			for(cItMap it=coeffsSmall_quasse.begin(); it != coeffsSmall_quasse.end(); ++it) {
				if(!isParallel && it->first.second > 1) continue; // If we are not in a parallel setting, discard parallel approximators

				// hard boundaries from observations
				if(isParallel && it->first.first == 1) continue;
				if(isParallel && nTips <= 64 && kStates <= 64 && it->first.second >= 2 ) continue;

				double predictedTime = predictLikelihoodEvaluationTimeQuasse(nTips, kStates, it->second);
				if(predictedTime < minTime) {
					minTime = predictedTime;
					bestSetting = it->first;
				}
			}
		} else { // large-scale model fitted
			for(cItMap it=coeffsLarge_quasse.begin(); it != coeffsLarge_quasse.end(); ++it) {
				if(!isParallel && it->first.second > 1) continue; // If we are not in a parallel setting, discard parallel approximators

				// hard boundaries from observations
				if(isParallel && it->first.first == 1) continue;
				if(isParallel && nTips <= 64 && kStates <= 64 && it->first.second >= 2 ) continue;

				double predictedTime = predictLikelihoodEvaluationTimeQuasse(nTips, kStates, it->second);
				if(predictedTime < minTime) {
					minTime = predictedTime;
					bestSetting = it->first;
				}
			}
		}
	} else if(ptrTensorCont->getOmega()->getDimensions().front() == 0) { // no cladogenetic events
		typedef std::map< std::pair<size_t, size_t>, std::vector<double> >::const_iterator cItMap;
		if(nTips <= thresholdNTips && kStates <= thresholdNStates) { // small-scale model fitted
			for(cItMap it=coeffsSmall.begin(); it != coeffsSmall.end(); ++it) {
				if(!isParallel && it->first.second > 1) continue; // If we are not in a parallel setting, discard parallel approximators

				// hard boundaries from observations
				if(it->first.second == 1 &&
				   (((it->first.first == 1 || it->first.first == 3) && nTips >= 64) ||
					((it->first.first == 2 || it->first.first == 4) && nTips <= 32))) {
					continue;
				}

				double predictedTime = predictLikelihoodEvaluationTimeMusse(nTips, kStates, it->second);
				if(predictedTime < minTime) {
					minTime = predictedTime;
					bestSetting = it->first;
				}
			}
		} else { // large-scale model fitted
			for(cItMap it=coeffsLarge.begin(); it != coeffsLarge.end(); ++it) {
				if(!isParallel && it->first.second > 1) continue; // If we are not in a parallel setting, discard parallel approximators

				// hard boundaries from observations
				if(it->first.second == 1 &&
				   (((it->first.first == 1 || it->first.first == 3) && nTips >= 64) ||
					((it->first.first == 2 || it->first.first == 4) && nTips <= 32))) {
					continue;
				}

				double predictedTime = predictLikelihoodEvaluationTimeMusse(nTips, kStates, it->second);
				if(predictedTime < minTime) {
					minTime = predictedTime;
					bestSetting = it->first;
				}
			}
		}
	} else if(ptrTensorCont->getOmega()->getDimensions()[0] > 0) { // cladogenetic events
		typedef std::map< std::pair<size_t, size_t>, std::vector<double> >::const_iterator cItMap;
		if(nTips <= thresholdNTips_clado && kStates <= thresholdNStates_clado) { // small-scale model fitted
			for(cItMap it=coeffsSmall_clado.begin(); it != coeffsSmall_clado.end(); ++it) {
				if(!isParallel && it->first.second > 1) continue; // If we are not in a parallel setting, discard parallel approximators

				// hard boundaries from observations
				if((it->first.first == 2 || it->first.first == 4) && nTips <= 40) { // We only consider compact/optimized for small trees
					continue;
				}

				double predictedTime = predictLikelihoodEvaluationTimeClasse(nTips, kStates, it->second);
				if(predictedTime < minTime) {
					minTime = predictedTime;
					bestSetting = it->first;
				}
			}
		} else { // large-scale model fitted
			for(cItMap it=coeffsLarge_clado.begin(); it != coeffsLarge_clado.end(); ++it) {
				if(!isParallel && it->first.second > 1) continue; // If we are not in a parallel setting, discard parallel approximators

				// hard boundaries from observations
				if((it->first.first == 2 || it->first.first == 4) && nTips <= 40) { // We only consider compact/optimized for small trees
					continue;
				}

				double predictedTime = predictLikelihoodEvaluationTimeClasse(nTips, kStates, it->second);
				if(predictedTime < minTime) {
					minTime = predictedTime;
					bestSetting = it->first;
				}
			}
		}
		if(isParallel && bestSetting.second > 1) {
			bestSetting.second = Utils::Parallel::Manager::getInstance()->getMaxNThread();
		}
	} else {
		assert(false && "Unknown case.");
	}

	// set nThreads
	size_t nThreads = std::min(bestSetting.second, Utils::Parallel::Manager::getInstance()->getMaxNThread());
	Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
	if(isParallel) {
		nThreadPerf.assign(Utils::Parallel::Manager::getInstance()->getMaxNThread()+1, accumulator_t(boost::accumulators::tag::rolling_window::window_size = ROLLING_MEAN_WINDOW_SIZE));
	}

	// return approximator
	return static_cast<approximatorVersion_t>(bestSetting.first);
}

ApproximatorSharedPtr AutoTuningApproximator::createApproximator(approximatorVersion_t anApproximatorType) {

	//std::cout << "Creating : " << (int)anApproximatorType << std::endl;
	if(anApproximatorType == SEQUENTIAL_OPTIMIZED) {
		return Factory::createSequentialTemplateCPU(intScheme, conditionType, ptrData, ptrScheduler, ptrSyncEventsCont, ptrTensorCont);
	} else if(anApproximatorType == SEQUENTIAL_BRANCHWISE) {
		return Factory::createSequentialBranchwiseCPU(intScheme, conditionType, ptrData, ptrScheduler, ptrSyncEventsCont, ptrTensorCont);
#if defined(_OPENMP)
	} else if(anApproximatorType == PARALLEL_OPTIMIZED) {
		return Factory::createParallelOpenMPCPU(intScheme, conditionType, ptrData, ptrScheduler, ptrSyncEventsCont, ptrTensorCont);
	} else if(anApproximatorType == PARALLEL_BRANCHWISE) {
		return Factory::createParallelBranchwiseCPU(intScheme, conditionType, ptrData, ptrScheduler, ptrSyncEventsCont, ptrTensorCont);
#endif
	} else {
		assert(false && "Unregistered approximator");
	}

	return ApproximatorSharedPtr();
}

} /* namespace Approximator */
} /* namespace Likelihood */
