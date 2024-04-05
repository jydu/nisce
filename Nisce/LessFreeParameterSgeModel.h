// SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * =====================================================================================
 *
 *       Filename:  LessFreeParameterSgeModel.h
 *
 * =====================================================================================
 */

#ifndef _LESS_FREE_PARAMETER_SGE_MODEL_H_
#define _LESS_FREE_PARAMETER_SGE_MODEL_H_

#include "SgeModel.h"

/**
 * @brief Free parameter model (mean and noise), but drop-out parameter has a priori distribution and is fitted with an empricial Bayesian approach.
 */
class LessFreeParameterSgeModel:
  public SgeModel2
{
public:
  string model_;

  // Distribution parameters
  size_t nbDistParams_;
  vector<string> distParamNames_;
  vector<shared_ptr<ConstraintInterface>> distParamConstraints_;
  d2_array_dbl distParamValues_;
  d3_array_dbl distArray_;    // Intermediate calculations. Contains, for each gene, each category, two numbers: value and probability.
  std::string prior_;

public:
  LessFreeParameterSgeModel(
      const d2_array_ui& counts,
      const d1_array_str& geneIds,
      const d1_array_str& cellIds,
      const string& model,
      size_t nbCategories,
      size_t nbRCategories,
      const std::string& prior = "Uniform") :
    SgeModel2(counts, geneIds, cellIds, nbCategories, nbRCategories, 2),
    model_(model),
    nbDistParams_(0),  // Will be initialized in the constructor body
    distParamNames_(), // Will be initialized in the constructor body
    distParamConstraints_(), // Will be initialized in the constructor body
    distParamValues_(nbGenes_),
    distArray_(nbGenes_),
    prior_(prior)
  {
    if (model == "ScaledBeta")
    {
      nbDistParams_ = 3;
      distParamNames_.resize(nbDistParams_);
      distParamNames_[0] = "Scale";
      distParamNames_[1] = "Alpha";
      distParamNames_[2] = "Beta";
      distParamConstraints_.resize(nbDistParams_);
      distParamConstraints_[0] = make_shared<IntervalConstraint>(0.01, 999, true, true);
      distParamConstraints_[1] = make_shared<IntervalConstraint>(0.01, 999, true, true);
      distParamConstraints_[2] = make_shared<IntervalConstraint>(0.01, 999, true, true);
    }
    else if (model == "Gamma")
    {
      nbDistParams_ = 2;
      distParamNames_.resize(nbDistParams_);
      distParamNames_[0] = "Mean";
      distParamNames_[1] = "Shape";
      distParamConstraints_.resize(nbDistParams_);
      distParamConstraints_[0] = make_shared<IntervalConstraint>(0.01, 999, true, true);
      distParamConstraints_[1] = make_shared<IntervalConstraint>(0.001, 999, true, true);
    }
    else
    {
      throw Exception("Unsupported model: " + model_);
    }
    for (auto& a : distParamValues_)
    {
      a.resize(nbDistParams_);
    }
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
  void computeCategories(size_t geneIndex, const d1_array_dbl& params)
  {
    if (model_ == "ScaledBeta")
    {
      computeCategoriesScaledBeta_(geneIndex, params);
      return;
    }
    if (model_ == "Gamma")
    {
      computeCategoriesGamma_(geneIndex, params);
      return;
    }
    throw Exception("Unsupported model: " + model_);
  }

  /**
   * @brief Set parameter values based on counts
   */
  void initParameters()
  {
    // First distribution parameter is the mean
    for (size_t i  = 0; i < nbGenes_; ++i)
    {
      size_t c = static_cast<size_t>(count(counts_[i].begin(), counts_[i].end(), 0));
      //Pseudo counts to avoid r == 0 or r == 1
      if (c == 0) c = 1;
      if (c == nbCells_) c = nbCells_ - 1;
      double s = accumulate(counts_[i].begin(), counts_[i].end(), 0.);
      distParamValues_[i][0] = s / static_cast<double>(nbCells_ - c);
      for (size_t j = 1; j < nbDistParams_; ++j)
      {
        distParamValues_[i][j] = 1.;
      }
    }
    // Library scale is 1.
    for (size_t j  = 0; j < nbCells_; ++j)
    {
      ls_[j] = 1.;
    }
    // Drop-out beta distribution parameters:
    if (prior_ == "Beta")
      computeRCategoriesBeta_({1., 1.});
    else if (prior_ == "Uniform")
      computeRCategoriesUniform_();
    else 
      throw Exception("Unknown prior distribution: " + prior_);
  }

  /**
   * @brief Output information about the model
   */
  void print(ostream& out) const
  {
    out << "-----------------------------------------------------------------" << endl;
    out << " SGE model with " << nbGenes_ << " genes and " << nbCells_ << " cells." << endl;
    out << " Poisson distribution is modulated by a " << model_ << " distribution," << endl;
    out << "     with " << nbCategories_ << " categories." << endl;
    out << " Drop-out parameters follow a " << prior_ << " prior distribution," << endl;
    out << "     with " << nbRCategories_ << " categories." << endl;
    out << "-----------------------------------------------------------------" << endl;
  }

  void computeRCategoriesBeta_(const d1_array_dbl& rHyperParams)
  {
    rHyperParams_ = rHyperParams;
    for (size_t i = 0; i < nbRCategories_; ++i)
    {
      rCatVals_[i] = boost::math::ibeta_inv(rHyperParams_[0], rHyperParams_[1], (static_cast<double>(i) + 0.5) / static_cast<double>(nbRCategories_)); // We use the median as category value
      rCatProbs_[i] = 1. / static_cast<double>(nbRCategories_);
    }
  }

  void computeRCategoriesUniform_()
  {
    for (size_t i = 0; i < nbRCategories_; ++i)
    {
      rCatVals_[i] = (static_cast<double>(i) + 0.5) / static_cast<double>(nbRCategories_);
      rCatProbs_[i] = 1. / static_cast<double>(nbRCategories_);
    }
  }

  /**
   * @brief Likelihood calculations
   */
  double computeLogLikelihoodForGene(size_t i, const d1_array_dbl& distParams)
  {
    computeCategories(i, distParams);
    double logLikelihood = 0.;
    for (size_t j = 0; j < nbCells_; ++j)
    {
      d1_array_dbl s(nbCategories_);
      for (size_t k = 0; k < nbCategories_; ++k)
      {
        s[k] = dpoisson(counts_[i][j], distArray_[i][k][0] * ls_[j]) * distArray_[i][k][1];
      }
      auto S = accumulate(s.begin(), s.end(), 0.);
      double P = 0;
      if (counts_[i][j] == 0)
      {
        for (size_t k = 0; k < nbRCategories_; ++k) {
          P += (rCatVals_[k] + (1. - rCatVals_[k]) * S) * rCatProbs_[k];
        }
      }
      else
      {
        for (size_t k = 0; k < nbRCategories_; ++k) {
          P += (1. - rCatVals_[k]) * S * rCatProbs_[k];
        }
      }
      logLikelihood += log(P);
    }
    if (isnan(logLikelihood) || isinf(logLikelihood))
    {
      throw (Exception("Error when computing likelihood."));
    }
    return logLikelihood;
  }

  double computeLogLikelihoodForCell(size_t j, double l)
  {
    double logLikelihood = 0.;
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      d1_array_dbl s(nbCategories_);
      for (size_t k = 0; k < nbCategories_; ++k)
      {
        s[k] = dpoisson(counts_[i][j], distArray_[i][k][0] * l) * distArray_[i][k][1];
      }
      auto S = accumulate(s.begin(), s.end(), 0.);
      double P = 0;
      if (counts_[i][j] == 0)
      {
        for (size_t k = 0; k < nbRCategories_; ++k) {
          P += (rCatVals_[k] + (1. - rCatVals_[k]) * S) * rCatProbs_[k];
        }
      }
      else
      {
        for (size_t k = 0; k < nbRCategories_; ++k) {
          P += (1. - rCatVals_[k]) * S * rCatProbs_[k];
        }
      }
      logLikelihood += log(P);
    }
    if (isnan(logLikelihood) || isinf(logLikelihood))
    { 
      throw (Exception("Error when computing likelihood."));
    }
    return logLikelihood;
  }

  double getLogLikelihoodForGene(size_t i) const override
  {
    return const_cast<LessFreeParameterSgeModel*>(this)->computeLogLikelihoodForGene(i, distParamValues_[i]); // Actually ok to do so here, since we are not modifying the object.
  }

  double getLogLikelihoodForCell(size_t j) const override
  {
    return const_cast<LessFreeParameterSgeModel*>(this)->computeLogLikelihoodForCell(j, ls_[j]); // Actually ok to do so here, since we are not modifying the object.
  }

  unsigned int getNumberOfParameters() const override
  {
    return static_cast<unsigned int>(nbCells_ + nbGenes_ * (nbDistParams_ + 1));
  }

  double getAicForGene(size_t i) const
  {
    return 2. * static_cast<double>(nbDistParams_ + 1) - 2. * getLogLikelihoodForGene(i); // Needs to account for r too (hence the +1).
  }

  double getMeanForGene(size_t i) const
  {
    if (model_ == "Gamma")
    {
      // This is the mean parameter:
      return distParamValues_[i][0];
    }
    else if (model_ == "ScaledBeta")
    {
      // This is scale * alpha / (alpha + beta)
      return distParamValues_[i][0] * distParamValues_[i][1] / (distParamValues_[i][1] + distParamValues_[i][2]);
    }
    throw Exception("Unsupported model: " + model_);
  }

  double getNoiseForGene(size_t i) const
  {
    if (model_ == "Gamma")
    {
      // This is 1/ shape
      return 1. / distParamValues_[i][1];
    }
    else if (model_ == "ScaledBeta")
    {
      // This is beta / (alpha * (alpha + beta + 1))
      return distParamValues_[i][2] / (distParamValues_[i][1] * (distParamValues_[i][1] + distParamValues_[i][2] + 1.));
    }
    throw Exception("Unsupported model: " + model_);
  }

  void printGeneParameters(ostream& out) const
  {
    out << "Gene,";
    for (auto p : distParamNames_)
    {
      out << "," << p;
    }
    out << ",logLk,nbParam,AIC" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i];
      for (auto p : distParamValues_[i])
      {
        out << "," << p;
      }
      out << "," << getLogLikelihoodForGene(i) << "," << static_cast<double>(nbDistParams_ + 1) << "," << getAicForGene(i) << endl;
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

  void printHyperParameters(ostream& out) const
  {
    out << "- Drop-out hyper parameters:" << endl;
    out << "  * Alpha: " << rHyperParams_[0] << endl;
    out << "  * Beta: "  << rHyperParams_[1] << endl;
  }

  void printResultsForGenes(ostream& out) const
  {
    out << "Gene,Mean,Noise" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i] << "," << getMeanForGene(i) << "," << getNoiseForGene(i) << endl;
    }
  }

private:
  void computeCategoriesScaledBeta_(size_t geneIndex, const d1_array_dbl& params)
  {
    for (size_t i = 0; i < nbCategories_; ++i)
    {
      distArray_[geneIndex][i][0] = params[0] * boost::math::ibeta_inv(params[1], params[2], ((double)i + 0.5) / (double)nbCategories_); // We use the median as category value
      distArray_[geneIndex][i][1] = 1. / static_cast<double>(nbCategories_);
    }
  }

  void computeCategoriesGamma_(size_t geneIndex, const d1_array_dbl& params)
  {
    for (size_t i = 0; i < nbCategories_; ++i)
    {
      distArray_[geneIndex][i][0] = params[0] * (1. / params[1]) * boost::math::gamma_p_inv(params[1], ((double)i + 0.5) / (double)nbCategories_); // We use the median as category value
      distArray_[geneIndex][i][1] = 1. / static_cast<double>(nbCategories_);
    }
  }

public:
  class GeneLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<LessFreeParameterSgeModel> sgeModel_;
    size_t geneId_;
    d1_array_dbl distParams_;

public:
    GeneLikelihood(shared_ptr<LessFreeParameterSgeModel> sgeModel, size_t geneId) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      geneId_(geneId),
      distParams_(sgeModel->nbDistParams_)
    {
      for (size_t i = 0; i < sgeModel->nbDistParams_; ++i)
      {
        addParameter_(new Parameter(
                        sgeModel->distParamNames_[i],
                        sgeModel->distParamValues_[geneId][i],
                        sgeModel->distParamConstraints_[i])
                      );
      }
    }

    GeneLikelihood* clone() const override { return new GeneLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      for (size_t j = 0; j < sgeModel_->nbDistParams_; ++j)
      {
        distParams_[j] = parameter(sgeModel_->distParamNames_[j]).getValue();
      }
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihoodForGene(geneId_, distParams_);
    }
  };

  class CellLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<LessFreeParameterSgeModel> sgeModel_;
    size_t cellId_;
    double l_;

public:
    CellLikelihood(shared_ptr<LessFreeParameterSgeModel> sgeModel, size_t cellId) :
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

  /**
   * @brief Used to estimate hyper-parameters
   */
  class FullLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<LessFreeParameterSgeModel> sgeModel_;
    d1_array_dbl params_;

public:
    FullLikelihood(shared_ptr<LessFreeParameterSgeModel> sgeModel) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      params_(sgeModel->rHyperParams_)
    {
      addParameter_(new Parameter("alpha", params_[0], Parameter::R_PLUS_STAR));
      addParameter_(new Parameter("beta", params_[1], Parameter::R_PLUS_STAR));
    }

    FullLikelihood* clone() const override { return new FullLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      params_[0] = parameter("alpha").getValue();
      params_[1] = parameter("beta").getValue();
      sgeModel_->computeRCategoriesBeta_(params_);
    }

    double getValue() const override
    {
      return -sgeModel_->getLogLikelihood();
    }
  };


};

#endif // _LESS_FREE_PARAMETER_SGE_MODEL_H_

