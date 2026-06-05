// SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * =====================================================================================
 *
 *       Filename:  EmpiricalBayesianSgeModel.h
 *
 * =====================================================================================
 */

#ifndef _EMPIRICAL_BAYESIAN_SGE_MODEL_H_
#define _EMPIRICAL_BAYESIAN_SGE_MODEL_H_


#include "SgeModel.h"

#include <limits>

// Only supports Gamma model. We assume all counts to share independent Gammas, but we could inforce a single one (like ASRV vs SSRV phylogenetic models).
class EmpiricalBayesianSgeModel:
  public SgeModel
{
public:
  size_t nbShapeCategories_;

  // Distribution parameters
  double meanShape_;
  double shapeShape_;

  d1_array_dbl distShapeValues_;

  d1_array_dbl means_;        // Gene-specific expression means
  d3_array_dbl distArray_;    // Intermediate calculations. Contains, for each shape category, each lambda category, two numbers: value and probability.

public:
  EmpiricalBayesianSgeModel(
      const d2_array_ui& counts,
      const d1_array_str& geneIds,
      const d1_array_str& cellIds,
      size_t nbShapeCategories,
      size_t nbCategories) :
    SgeModel(counts, geneIds, cellIds, nbCategories),
    nbShapeCategories_(nbShapeCategories),
    meanShape_(1.),
    shapeShape_(1.),
    distShapeValues_(nbShapeCategories),
    means_(nbGenes_),
    distArray_(nbShapeCategories_)
  {
    for (auto& a : distArray_)
    {
      a.resize(nbCategories_);
      for (auto& b : a)
      {
        b.resize(2);
      }
    }
  }

public:
  /**
   * @brief Recomputes the Discretized Poisson modulated PDF for a given gene.
   */
  void computeCategories()
  {
    for (size_t j = 0; j < nbShapeCategories_; ++j)
    {
      distShapeValues_[j] = meanShape_ * (1. / shapeShape_) * boost::math::gamma_p_inv(shapeShape_, ((double)j + 0.5) / (double)nbCategories_); // We use the median as category value
    }
    for (size_t j = 0; j < nbShapeCategories_; ++j)
    {
      for (size_t k = 0; k < nbCategories_; ++k)
      {
        distArray_[j][k][0] = (1. / distShapeValues_[j]) * boost::math::gamma_p_inv(distShapeValues_[j], ((double)k + 0.5) / (double)nbCategories_); // We use the median as category value
        distArray_[j][k][1] = 1. / static_cast<double>(nbCategories_);
      }
    }
  }

  /**
   * @brief Set parameter values based on counts
   */
  void initParameters()
  {
    // Dropout probability is proportion of 0s
    for (size_t i  = 0; i < nbGenes_; ++i)
    {
      size_t c = static_cast<size_t>(count(counts_[i].begin(), counts_[i].end(), 0));
      //Pseudo counts to avoid r == 0 or r == 1
      if (c == 0) c = 1;
      if (c == nbCells_) c = nbCells_ - 1;
      double s = accumulate(counts_[i].begin(), counts_[i].end(), 0.);
      rs_[i] = static_cast<double>(c) / static_cast<double>(nbCells_);
      means_[i] = s / static_cast<double>(nbCells_ - c);
    }
    // Library scale is 1.
    for (size_t j  = 0; j < nbCells_; ++j)
    {
      ls_[j] = 1.;
    }
  }

  /**
   * @brief Output information about the model
   */
  void print(ostream& out) const
  {
    out << "-----------------------------------------------------------------" << endl;
    out << " SGE model with " << nbGenes_ << " genes and " << nbCells_ << " cells." << endl;
    out << " Poisson distribution is modulated by a gamma distribution." << endl;
    out << " Noise param. are estimated via an empirical Bayesian approach." << endl;
    out << "-----------------------------------------------------------------" << endl;
  }

  /**
   * @brief Likelihood calculations
   */
  double computeLogLikelihoodForGeneForClass(size_t i, double r, double mean, size_t shapeClass)
  {
    double conditionalLogLikelihood = 0.;
    for (size_t j = 0; j < nbCells_; ++j)
    {
      double l = ls_[j];
      d1_array_dbl s(nbCategories_);
      for (size_t k = 0; k < nbCategories_; ++k)
      {
        s[k] = dpoisson(counts_[i][j], mean * l * distArray_[shapeClass][k][0]) * distArray_[shapeClass][k][1];
      }
      auto S = accumulate(s.begin(), s.end(), 0.);
      double P;
      if (counts_[i][j] == 0)
      {
        P = r + (1. - r) * S;
      }
      else
      {
        P = (1. - r) * S;
      }
      if (isnan(log(P)) || isinf(log(P)))
      {
        P = numeric_limits<double>::min();
      }
      conditionalLogLikelihood += log(P);
    }
    if (isnan(conditionalLogLikelihood) || isinf(conditionalLogLikelihood))
    {
      throw (Exception("Conditional likelihood is 0!!!"));
    }
    return conditionalLogLikelihood;
  }

  double computeLogLikelihoodForGene(size_t i, double r, double mean)
  {
    d1_array_dbl logLikelihoods(nbShapeCategories_);
    for (size_t k = 0; k < nbShapeCategories_; ++k)
    {
      logLikelihoods[k] = computeLogLikelihoodForGeneForClass(i, r, mean, k) - log(static_cast<double>(nbShapeCategories_));
    }
    double logLikelihood = logSumExp(logLikelihoods.begin(), logLikelihoods.end());
    if (isnan(logLikelihood) || isinf(logLikelihood))
    {
      throw (Exception("Likelihood is 0!!!"));
    }
    return logLikelihood;
  }

  //an alternative version where genes have independent ditributions per gene
  double computeLogLikelihoodForGene2(size_t i, double r, double mean)
  {
    double logLikelihood = 0.;
    for (size_t j = 0; j < nbCells_; ++j)
    {
      double l = ls_[j];
      d1_array_dbl s(nbShapeCategories_ * nbCategories_);
      size_t z = 0;
      for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
      {
        for (size_t k3 = 0; k3 < nbCategories_; ++k3)
        {
          s[z++] = dpoisson(counts_[i][j], mean * l * distArray_[k2][k3][0]) * distArray_[k2][k3][1];
        }
      }
      auto S = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbShapeCategories_);
      double P;
      if (counts_[i][j] == 0)
      {
        P = r + (1. - r) * S;
      }
      else
      {
        P = (1. - r) * S;
      }
      logLikelihood += log(P);
    }
    if (isnan(logLikelihood) || isinf(logLikelihood))
    {
      throw (Exception("Likelihood is 0!!!"));
    }
    return logLikelihood;
  }

  double computeLogLikelihoodForCell(size_t j, double l)
  {
    double logLikelihood = 0.;
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      double mean = means_[i];
      d1_array_dbl s(nbShapeCategories_ * nbCategories_);
      size_t z = 0;
      for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
      {
        for (size_t k3 = 0; k3 < nbCategories_; ++k3)
        {
          s[z++] = dpoisson(counts_[i][j], mean * l * distArray_[k2][k3][0]) * distArray_[k2][k3][1];
        }
      }
      auto S = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbShapeCategories_);
      double P;
      double r = rs_[i];
      if (counts_[i][j] == 0)
      {
        P = r + (1. - r) * S;
      }
      else
      {
        P = (1. - r) * S;
      }
      logLikelihood += log(P);
    }
    if (isnan(logLikelihood) || isinf(logLikelihood))
    {
      throw (Exception("Likelihood is 0!!!"));
    }
    return logLikelihood;
  }

  double computeLogLikelihood(double meanShape, double shapeShape)
  {
    meanShape_ = meanShape;
    shapeShape_ = shapeShape;
    computeCategories();
    d1_array_dbl logLik(nbGenes_);

    for_each (std::execution::par, geneIndex_.begin(), geneIndex_.end(), [&](size_t i) {
      logLik[i] = computeLogLikelihoodForGene(i, rs_[i], means_[i]);
    });
    sort(logLik.begin(), logLik.end()); // Sort by increasing value to limit underflow.
    double logLikelihood = accumulate(logLik.begin(), logLik.end(), 0.);

    return logLikelihood;
  }

  unsigned int getNumberOfParametersForGene(size_t i) const override
  {
    return 2;
  }

  unsigned int getNumberOfParameters() const override
  {
    return static_cast<unsigned int>(nbCells_ + nbGenes_ * 2 + 2); //1 l per cell, 1 r and 1 mean per gene, 2 hyper parameters for noise.
  }

  double getLogLikelihoodForGene(size_t i) const override
  {
    return const_cast<EmpiricalBayesianSgeModel*>(this)->computeLogLikelihoodForGene(i, rs_[i], means_[i]); // Actually ok to do so here, since we are not modifying the object.
  }

  double getLogLikelihoodForCell(size_t j) const override
  {
    return const_cast<EmpiricalBayesianSgeModel*>(this)->computeLogLikelihoodForCell(j, ls_[j]); // Actually ok to do so here, since we are not modifying the object.
  }

  double getMeanForGene(size_t i) const
  {
    return means_[i];
  }

  NumberResult<double> getNoiseForGene(size_t i) const
  {
    d1_array_dbl logShapeLikelihoods(nbShapeCategories_);
    for (size_t k = 0; k < nbShapeCategories_; ++k)
    {
      logShapeLikelihoods[k] = const_cast<EmpiricalBayesianSgeModel*>(this)->computeLogLikelihoodForGeneForClass(i, rs_[i], means_[i], k) - log(static_cast<double>(nbShapeCategories_));
    }
    double logLik = logSumExp(logShapeLikelihoods.begin(), logShapeLikelihoods.end());
    double postNoise = 0.;
    double entropyNoise = 0.;
    for (size_t k = 0; k < nbShapeCategories_; ++k)
    {
      double p = exp(logShapeLikelihoods[k] - logLik);
      if (p > 0) {
        postNoise += p * (1. / distShapeValues_[k]);
        entropyNoise -= p * log(p);
      }
    }
    return NumberResult<double>(postNoise, entropyNoise);
  }

 
  void printHyperParameters(ostream& out) const
  {
    out << "- Shape hyper parameters:" << endl;
    out << "  * Mean: " << meanShape_ << endl;
    out << "  * Shape: " << shapeShape_ << endl;
  }

  void printGeneParameters(ostream& out) const
  {
    out << "Gene,r,mean,logLk" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i] << "," << rs_[i] << "," << means_[i] << "," << getLogLikelihoodForGene(i) << endl;
    }
  }

  void printCellParameters(ostream& out) const
  {
    out << "Cell,l" << endl;
    for (size_t j = 0; j < nbCells_; ++j)
    {
      out << cellIds_[j] << "," << ls_[j] << endl;
    }
  }

  void printResultsForGenes(ostream& out) const
  {
    out << "Gene,Mean,Noise,EntropyNoise" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      double resMean = getMeanForGene(i);
      auto resNoise = getNoiseForGene(i);
      out << geneIds_[i] << "," << resMean << "," << resNoise.value << "," << resNoise.score << endl;
    }
  }


public:
  class EBLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<EmpiricalBayesianSgeModel> sgeModel_;
    double meanShape_;
    double shapeShape_;

public:
    EBLikelihood(shared_ptr<EmpiricalBayesianSgeModel> sgeModel) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      meanShape_(0),
      shapeShape_(0)
    {
      addParameter_(new Parameter("meanShape", sgeModel->meanShape_, make_shared<IntervalConstraint>(0.01, 999, true, true)));
      addParameter_(new Parameter("shapeShape", sgeModel->shapeShape_, make_shared<IntervalConstraint>(0.01, 999, true, true)));
    }

    EBLikelihood* clone() const override { return new EBLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      meanShape_ = parameter("meanShape").getValue();
      shapeShape_ = parameter("shapeShape").getValue();
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihood(meanShape_, shapeShape_);
    }
  };

  class GeneLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<EmpiricalBayesianSgeModel> sgeModel_;
    size_t geneId_;
    double r_;
    double mean_;

public:
    GeneLikelihood(shared_ptr<EmpiricalBayesianSgeModel> sgeModel, size_t geneId) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      geneId_(geneId),
      r_(0),
      mean_(1.)
    {
      addParameter_(new Parameter("r", sgeModel->rs_[geneId], Parameter::PROP_CONSTRAINT_EX));
      addParameter_(new Parameter("mean", sgeModel->rs_[geneId], Parameter::R_PLUS_STAR));
    }

    GeneLikelihood* clone() const override { return new GeneLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      r_ = parameter("r").getValue();
      mean_ = parameter("mean").getValue();
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihoodForGene(geneId_, r_, mean_);
    }
  };

  class CellLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<EmpiricalBayesianSgeModel> sgeModel_;
    size_t cellId_;
    double l_;

public:
    CellLikelihood(shared_ptr<EmpiricalBayesianSgeModel> sgeModel, size_t cellId) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      cellId_(cellId),
      l_(0)
    {
      addParameter_(new Parameter("l", sgeModel->ls_[cellId], Parameter::R_PLUS_STAR));
    }

    CellLikelihood* clone() const override { return new CellLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      l_ = parameter("l").getValue();
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihoodForCell(cellId_, l_);
    }
  };
};

#endif //_EMPIRICAL_BAYESIAN_SGE_MODEL_H_

