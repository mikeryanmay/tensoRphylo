/*
 * Utils.h
 *
 *  Created on: Sep 16, 2019
 *      Author: xaviermeyer
 */

#ifndef TEST_UTILS_H_
#define TEST_UTILS_H_

#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Tensor/IncTensor.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Parameters/Container.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Likelihood/Approximator/IncLikelihoodApproximator.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.h"
#include "Utils/MemoryPool/EigenCPU.h"
#include "Likelihood/ConditionTypes/ConditionType.h"
#include "Utils/Output/OutputManager.h"

#include <boost/assign/list_of.hpp>

#include "../Interface/DistributionHandler.h"
#include "../Likelihood/Approximator/StochasticMapping/IncFwdStochasticMapping.h"
namespace Test {

namespace Definitions {
	static const std::vector< std::string > intSchemeNames = boost::assign::list_of("EULER")("RUNGE_KUTTA4")("RUNGE_KUTTA54")("RUNGE_KUTTA_DOPRI5");
} /* namespace Definitions */

namespace InitData {

	Phylogeny::NexusReader::NexusParserSharedPtr nexusFromFile(std::string fileName);
	Phylogeny::NexusReader::NexusParserSharedPtr nexusFromFile(std::string nexusFilename, std::string tsvFilename);
	Phylogeny::NexusReader::NexusParserSharedPtr nexusSingleBranch(double branchLength, size_t nState, std::vector<int> &myStates);
	Phylogeny::NexusReader::NexusParserSharedPtr nexusFromSingleBranchFile(std::string fileName);

	Parameters::ContainerSharedPtr parametersFromFile(std::string fileName);
	Parameters::ContainerSharedPtr parametersYule1(size_t nState, double lambda, double eta, double rho);
	Parameters::ContainerSharedPtr parametersBD1(size_t nState, double lambda, double mu, double eta, double rho);
	Parameters::ContainerSharedPtr parametersFBD1(size_t nState, double lambda, double mu, double phi, double eta, double rho);
	Parameters::ContainerSharedPtr parametersBDME(size_t nState, double lambda, double mu, double eta, double rho, double me_time, double me_prob);
	Parameters::ContainerSharedPtr parametersBiSSE(size_t nState, std::vector<double> lambda, std::vector<double> mu, double eta, double rho);


} /* namespace InitData */

namespace SequentialCPU {
	double computeLikForTests(Likelihood::Integrator::integrationScheme_t intScheme,
					  Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					  Parameters::ContainerSharedPtr parameters);

	double computeLikForInterface(Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
								  const std::string &parametersFile);

	TensorPhylo::Interface::vecHistories_t computeStochasticMappingForInterface(
			size_t nHistories,
			Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
			const std::string &parametersFile);

	double computeLikForInterfaceTests(int idApproximator,
									   Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
									   const std::string &parametersFile);

	double computeLikForInterfaceTests(Likelihood::Integrator::integrationScheme_t intScheme,
			   	   	   const std::string &aNewickStringRB,
					   Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					   Parameters::ContainerSharedPtr parameters);

	double computeLogLik(size_t nLikApproximation,
						 Likelihood::Integrator::integrationScheme_t intScheme,
						 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
						 Parameters::ContainerSharedPtr parameters);

	TensorPhylo::Interface::vecHistories_t computeStochasticMappingLegacy(size_t nReplicas,
						 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
						 Parameters::ContainerSharedPtr parameters,
						 PS::TreeSharedPtr ptrTree);

	TensorPhylo::Interface::vecHistories_t computeStochasticMappingLegacy(size_t nReplicas,
	 	 	 	 	 	 Likelihood::Approximator::StochasticMapping::stochasticMappingAlgo_t aAlgo,
						 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
						 Parameters::ContainerSharedPtr parameters,
						 PS::TreeSharedPtr ptrTree);


	void runBenchmarkLegacy(size_t nLikApproximation,
						Likelihood::Integrator::integrationScheme_t intScheme,
						Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
						Parameters::ContainerSharedPtr parameters,
						const std::string &logFile);

	void runBenchmarkNew(size_t nLikApproximation,
						Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
						const std::string &parametersFile,
						const std::string &logFile);

	void benchmarkConvolution(const size_t M, const std::vector<size_t> &vecK);
	void benchmarkConvolution2(const size_t M, const std::vector<size_t> &vecK);

	bool checkHistoriesEndStatesAgainstDataProbability(
			PS::TreeSharedPtr ptrTree,
			Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
			TensorPhylo::Interface::vecHistories_t vecHistories,
			double tolerance);

} /* namespace SequentialCPU */

} /* namespace Test */

#endif /* TEST_UTILS_H_ */
