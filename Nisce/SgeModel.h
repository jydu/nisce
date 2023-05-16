/*
 * =====================================================================================
 *
 *       Filename:  SgeModel.h
 *
 * =====================================================================================
 */

#ifndef _SGE_MODEL_H_
#define _SGE_MODEL_H_

#include "ProbFunctions.h"

class SgeModel
{
public:
  d2_array_ui counts_;
  d1_array_str geneIds_;
  d1_array_str cellIds_;
  size_t nbGenes_;
  size_t nbCells_;
  size_t nbCategories_;

  d1_array_dbl rs_;           // Dropout parameters
  d1_array_dbl ls_;           // Library size parameters

  vector<size_t> geneIndex_;
  vector<size_t> cellIndex_;

public:
  SgeModel(
    const d2_array_ui& counts,
    const d1_array_str& geneIds,
    const d1_array_str& cellIds,
    size_t nbCategories) :
    counts_(counts),
    geneIds_(geneIds),
    cellIds_(cellIds),
    nbGenes_(counts.size()),
    nbCells_(nbGenes_ == 0 ? 0 : counts[0].size()),
    nbCategories_(nbCategories),
    rs_(nbGenes_),
    ls_(nbCells_),
    geneIndex_(),
    cellIndex_()
  {
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      geneIndex_.push_back(i);
    }
    for (size_t j = 0; j < nbCells_; ++j)
    {
      cellIndex_.push_back(j);
    }
  }

  virtual ~SgeModel() {}

public:
  virtual unsigned int getNumberOfParameters() const = 0;

  virtual double getLogLikelihoodForGene(size_t i) const = 0;
  
  virtual double getLogLikelihoodForCell(size_t j) const = 0;

  virtual double getLogLikelihood() const
  {
    // We use the cell likelihood, has it does not recomputes the categories and is then faster
    d1_array_dbl logCellLikelihoods(nbCells_);
    for_each (std::execution::par,
	      cellIndex_.begin(),
	      cellIndex_.end(),
	      [&](size_t j)
    {
      logCellLikelihoods[j] = getLogLikelihoodForCell(j);
    });
    // sort to limit underflow:
    std::sort(logCellLikelihoods.begin(), logCellLikelihoods.end());
    double logLikelihood = std::accumulate(logCellLikelihoods.begin(), logCellLikelihoods.end(), 0.);
    return logLikelihood;
  }

};

#endif //_SGE_MODEL_H_

