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

/**
 * @brief Basal interface for SGE models
 */
class SgeModelInterface
{
public:
  virtual ~SgeModelInterface() {}
  
  virtual unsigned int getNumberOfParameters() const = 0;

  virtual double getLogLikelihoodForGene(size_t i) const = 0;
  
  virtual double getLogLikelihoodForCell(size_t j) const = 0;

  virtual double getLogLikelihood() const = 0;

};

/**
 * @brief Partial implementation of an SgeModel interface where drop-out and normalization parameters are free.
 */
class SgeModel:
  public virtual SgeModelInterface
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

  virtual double getLogLikelihood() const override
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


/**
 * @brief Partial implementation of the SgeModel interface when drop-out parameters have an a priori distribution.
 */
class SgeModel2:
  public virtual SgeModelInterface
{
public:
  d2_array_ui counts_;
  d1_array_str geneIds_;
  d1_array_str cellIds_;
  size_t nbGenes_;
  size_t nbCells_;
  size_t nbCategories_;
  size_t nbRCategories_;
  size_t nbRHyperParams_;

  d1_array_dbl rCatVals_;     // Dropout cat values
  d1_array_dbl rCatProbs_;    // Dropout cat probabilities
  d1_array_dbl ls_;           // Library size parameters
  d1_array_dbl rHyperParams_;
  
  vector<size_t> geneIndex_;
  vector<size_t> cellIndex_;

public:
  SgeModel2(
    const d2_array_ui& counts,
    const d1_array_str& geneIds,
    const d1_array_str& cellIds,
    size_t nbCategories,
    size_t nbRCategories,
    size_t nbRHyperParams) :
    counts_(counts),
    geneIds_(geneIds),
    cellIds_(cellIds),
    nbGenes_(counts.size()),
    nbCells_(nbGenes_ == 0 ? 0 : counts[0].size()),
    nbCategories_(nbCategories),
    nbRCategories_(nbRCategories),
    nbRHyperParams_(nbRHyperParams),
    rCatVals_(nbRCategories_),
    rCatProbs_(nbRCategories_),
    ls_(nbCells_),
    rHyperParams_(nbRHyperParams),
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

  virtual ~SgeModel2() {}

public:
  virtual double getLogLikelihood() const override
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

