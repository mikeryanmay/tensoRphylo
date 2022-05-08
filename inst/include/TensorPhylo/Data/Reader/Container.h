/*
 * Container.h
 *
 *  Created on: Mar 10, 2020
 *      Author: meyerx
 */

#ifndef DATA_READER_CONTAINER_H_
#define DATA_READER_CONTAINER_H_

#include <cstddef>
#include <Eigen/Core>

#include <map>
#include <string>
#include <vector>

namespace Phylogeny {
namespace Data {

class Container {
public:
	Container();
	~Container();

	bool isReady() const;

	void setTaxaToIdMap(const std::vector<std::string> &aTaxaVec);
	void setTaxaMaps(const std::map<std::string, std::size_t> &aTaxaToIdMap, const std::map<std::size_t, std::string> &aIdToTaxaMap);

	void registerTaxaData(const std::string &taxaName, const std::vector<double> &probs);
	void registerTaxaData(size_t nStates, const std::string &taxaName, const std::vector<int> &states);

	const Eigen::VectorXd& getProbForTaxaId(int idTaxa) const;
	const Eigen::VectorXd& getProbForTaxaIdThreadSafe(int idTaxa) const;
	const Eigen::VectorXd& getProbForTaxaLabel(const std::string &taxaLabel) const;

private:

	std::map<size_t, std::string> mapIdToTaxa;
	std::map<std::string, size_t> mapTaxaToId;
	std::map< std::string, Eigen::VectorXd  > mapTaxaToProbs;

	typedef struct  {
		int idTaxa;
		Eigen::VectorXd statesProb;
	} stateCache_t;

	mutable std::vector<stateCache_t> vecProbCache;

	size_t getPositionForTaxaId(int idTaxa) const;

};

} /* namespace Data */
} /* namespace Phylogeny */

#endif /* DATA_READER_CONTAINER_H_ */
