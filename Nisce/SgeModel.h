// SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
//
// SPDX-License-Identifier: GPL-3.0-or-later

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
  
  virtual unsigned int getNumberOfParametersForGene(size_t i) const = 0;
  
  virtual unsigned int getNumberOfParameters() const = 0;

  virtual double getLogLikelihoodForGene(size_t i) const = 0;
  
  virtual double getLogLikelihoodForCell(size_t j) const = 0;

  virtual double getLogLikelihood() const = 0;

  virtual size_t getNumberOfGenes() const = 0;
  
  virtual size_t getNumberOfCells() const = 0;

};

/**
 * @brief Basal interface for SGE models
 */
class SgeSimulatorInterface
{
public:
  virtual ~SgeSimulatorInterface() {}
  
  virtual unsigned int drawCount(size_t geneIndex, size_t cellIndex) const = 0;
  
  virtual std::unique_ptr<vector<vector<unsigned int>>> simulateCounts() const = 0;

  virtual size_t getNumberOfGenes() const = 0;
  
  virtual size_t getNumberOfCells() const = 0;

};

/**
 * @brief Basal interface for SGE models
 */
class AbstractSgeSimulator:
  public virtual SgeSimulatorInterface
{
public:
  virtual ~AbstractSgeSimulator() {}
  
  std::unique_ptr<d2_array_ui> simulateCounts() const override
  {
    auto counts = std::make_unique<d2_array_ui>(getNumberOfGenes());
    for (size_t i = 0; i < getNumberOfGenes(); ++i) {
      (*counts)[i].resize(getNumberOfCells());
      for (size_t j = 0; j < getNumberOfCells(); ++j) {
	(*counts)[i][j] = drawCount(i, j);
      }
    }
    return counts;
  }

};

/**
 * @brief Partial implementation of an SgeModel interface where drop-out and normalization parameters are free.
 */
class SgeModel:
  public virtual SgeModelInterface,
  public AbstractSgeSimulator
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
  bool withDropout_;

  vector<size_t> geneIndex_;
  vector<size_t> cellIndex_;

public:
  SgeModel(
    const d2_array_ui& counts,
    const d1_array_str& geneIds,
    const d1_array_str& cellIds,
    size_t nbCategories,
    bool withDropout) :
    counts_(counts),
    geneIds_(geneIds),
    cellIds_(cellIds),
    nbGenes_(counts.size()),
    nbCells_(nbGenes_ == 0 ? 0 : counts[0].size()),
    nbCategories_(nbCategories),
    rs_(nbGenes_),
    ls_(nbCells_),
    withDropout_(withDropout),
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
  
  size_t getNumberOfGenes() const override { return nbGenes_; }
  
  size_t getNumberOfCells() const override { return nbCells_; }

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

  unsigned int drawCount(size_t geneIndex, size_t cellIndex) const override { return 0; }
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
  size_t getNumberOfGenes() const override { return nbGenes_; }
  
  size_t getNumberOfCells() const override { return nbCells_; }

  virtual double getLogLikelihood() const override
  {
    // As r is shared by gene, we use the gene likelihood.
    d1_array_dbl logGeneLikelihoods(nbGenes_);
    for_each (std::execution::par,
	      geneIndex_.begin(),
	      geneIndex_.end(),
	      [&](size_t i)
    {
      logGeneLikelihoods[i] = getLogLikelihoodForGene(i);
    });
    // sort to limit underflow:
    std::sort(logGeneLikelihoods.begin(), logGeneLikelihoods.end());
    double logLikelihood = std::accumulate(logGeneLikelihoods.begin(), logGeneLikelihoods.end(), 0.);
    return logLikelihood;
  }

};

class OptimizationParameters
{
public:
  std::string description;
  double tolerance;
  std::string outputGenes;
  std::string outputCells;
  std::string outputGenesHyper;
  bool estimateNorms;
  bool estimateHyperParameters;
  bool aliasShape;
};
#endif //_SGE_MODEL_H_

