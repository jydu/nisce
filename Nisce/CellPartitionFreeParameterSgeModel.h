// SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * =====================================================================================
 *
 *       Filename:  CellPartitionFreeParameterSgeModel.h
 *
 * =====================================================================================
 */

#ifndef _CELL_PARTITION_FREE_PARAMETER_SGE_MODEL_H_
#define _CELL_PARTITION_FREE_PARAMETER_SGE_MODEL_H_

#include "SgeModel.h"
#include "Partition.h"

// From the STL:
#include <map>

// From bpp-core:
#include <Bpp/Numeric/AbstractParameterAliasable.h>

// In this model, there is one mean and noise parameter per partition.
class CellPartitionFreeParameterSgeModel:
  public SgeModel,
  public enable_shared_from_this<CellPartitionFreeParameterSgeModel>
{
public:
  Partition partition_;
  
  // Distribution parameters
  size_t nbDistParams_;
  std::vector<string> distParamNames_;
  std::vector<std::shared_ptr<ConstraintInterface>> distParamConstraints_;
  d2_array_dbl distParamValues_;
  d4_array_dbl distArray_;    // Intermediate calculations. Contains, for each gene, each partition elt, each category, two numbers: value and probability.
  bool shapeIsAliased_;       // Record whether shapes have been aliased during optimization to get the number of parameters right.

public:
  CellPartitionFreeParameterSgeModel(
      const d2_array_ui& counts,
      const d1_array_str& geneIds,
      const d1_array_str& cellIds,
      const Partition& partition,
      size_t nbCategories) :
    SgeModel(counts, geneIds, cellIds, nbCategories),
    partition_(partition),
    nbDistParams_(2 * partition.size()),
    distParamNames_(), // Will be initialized in the constructor body
    distParamConstraints_(), // Will be initialized in the constructor body
    distParamValues_(nbGenes_),
    distArray_(nbGenes_),
    shapeIsAliased_(false)
  {
    for (auto& p : partition_.partitionMap) {
      distParamNames_.push_back(p.first + ".Mean");
      distParamNames_.push_back(p.first + ".Shape");
      distParamConstraints_.push_back(make_shared<IntervalConstraint>(0.01, 999, true, true));
      distParamConstraints_.push_back(make_shared<IntervalConstraint>(0.01, 999, true, true));
    }
    for (auto& a : distParamValues_)
    {
      a.resize(nbDistParams_);
    }
    for (auto& a : distArray_)
    {
      a.resize(partition_.size());
      for (auto& b : a)
      {
        b.resize(nbCategories_);
        for (auto& c : b)
        {
          c.resize(2);
        }
      }
    }
  }

public:
  /**
   * @brief Recomputes the Discretized Poisson modulated PDF for a given gene.
   */
  void computeCategories(size_t geneIndex, size_t partitionIndex, double mean, double shape)
  {
    for (size_t k = 0; k < nbCategories_; ++k)
    {
      distArray_[geneIndex][partitionIndex][k][0] = mean * (1. / shape) * boost::math::gamma_p_inv(shape, ((double)k + 0.5) / (double)nbCategories_); // We use the median as category value
      distArray_[geneIndex][partitionIndex][k][1] = 1. / static_cast<double>(nbCategories_);
    }
  }

  /**
   * @brief Set parameter values based on counts
   */
  void initParameters()
  {
    // Dropout probability is proportion of 0s
    // First distribution parameter is the mean
    for (size_t i  = 0; i < nbGenes_; ++i)
    {
      size_t c = static_cast<size_t>(count(counts_[i].begin(), counts_[i].end(), 0));
      //Pseudo counts to avoid r == 0 or r == 1
      if (c == 0) c = 1;
      if (c == nbCells_) c = nbCells_ - 1;
      double s = accumulate(counts_[i].begin(), counts_[i].end(), 0.);
      rs_[i] = static_cast<double>(c) / static_cast<double>(nbCells_);
      for (size_t p = 0; p < partition_.size(); ++p)
      {
        distParamValues_[i][2 * p] = s / static_cast<double>(nbCells_ - c);
        distParamValues_[i][2 * p + 1] = 1.;
      }
    }
    // Library scale is the library size divided by the mean library size of the partition.
    for (const auto& p : partition_.partitionMap)
    {
      double meanLibSize = 0.;
      for (auto j : p.second)
      {
        unsigned int libSize = 0;
        for (size_t i  = 0; i < nbGenes_; ++i)
        {
          libSize += counts_[i][j];
        }
        ls_[j] = static_cast<double>(libSize);
        meanLibSize += static_cast<double>(libSize) / static_cast<double>(p.second.size());
      }
      for (auto j : p.second)
      {
        ls_[j] /= meanLibSize;
      }
    }
  }

  /**
   * @brief Output information about the model
   */
  void print(ostream& out) const
  {
    out << "-----------------------------------------------------------------" << endl;
    out << " SGE model with " << nbGenes_ << " genes and " << nbCells_ << " cells." << endl;
    out << " Model has " << partition_.size() << " cell groups." << endl;
    out << " Poisson distribution is modulated by a Gamma distribution," << endl;
    out << "     with " << nbCategories_ << " categories." << endl;
    out << "-----------------------------------------------------------------" << endl;
  }

  /**
   * @brief Likelihood calculations
   */
  double computeLogLikelihoodForGene(size_t i, double r, const d1_array_dbl& distParams)
  {
    for (size_t p = 0; p < partition_.size(); ++p)
    {
      computeCategories(i, p, distParams[2 * p], distParams[2 * p + 1]);
    }
    double logLikelihood = 0.;
    size_t partitionIndex = 0;
    for (auto& p : partition_.partitionMap)
    {
      for (size_t j : p.second)
      {
        d1_array_dbl s(nbCategories_);
        for (size_t k = 0; k < nbCategories_; ++k)
        {
          s[k] = dpoisson(counts_[i][j], distArray_[i][partitionIndex][k][0] * ls_[j]) * distArray_[i][partitionIndex][k][1];
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
        logLikelihood += log(P);
      }
      partitionIndex++;
    }
    if (isnan(logLikelihood) || isinf(logLikelihood))
    {
      throw (Exception("Conditional likelihood is 0!!!"));
    }
    return logLikelihood;
  }
  
  double computeLogLikelihoodForCell(size_t j, double l)
  {
    double logLikelihood = 0.;
    size_t partitionIndex = partition_.getPartitionIndex(j);
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      d1_array_dbl s(nbCategories_);
      for (size_t k = 0; k < nbCategories_; ++k)
      {
        s[k] = dpoisson(counts_[i][j], distArray_[i][partitionIndex][k][0] * l) * distArray_[i][partitionIndex][k][1];
      }
      auto S = accumulate(s.begin(), s.end(), 0.);
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
      throw (Exception("Conditional likelihood is 0!!!"));
    }
    return logLikelihood;
  }

  double getLogLikelihoodForGene(size_t i) const override
  {
    return const_cast<CellPartitionFreeParameterSgeModel*>(this)->computeLogLikelihoodForGene(i, rs_[i], distParamValues_[i]); // Actually ok to do so here, since we are not modifying the object.
  }

  double getLogLikelihoodForCell(size_t j) const override
  {
    return const_cast<CellPartitionFreeParameterSgeModel*>(this)->computeLogLikelihoodForCell(j, ls_[j]); // Actually ok to do so here, since we are not modifying the object.
  }

  unsigned int getNumberOfParametersForGene(size_t i) const override
  {
    //Note: all genes have the same number of parameters
    unsigned int np = 0;
    if(shapeIsAliased_) {
      np = static_cast<unsigned int>((nbDistParams_ - partition_.size()) + 2);
      //                             |      mean only                  |   |
      //                                                                   | r and shape
      //
    } else {
      np =  static_cast<unsigned int>(nbDistParams_ + 1);
      //                            |mean and shape| |
      //                                             |r
    }
    return np;    
  }

  unsigned int getNumberOfParameters() const override
  {
    return static_cast<unsigned int>(nbCells_ + nbGenes_ * getNumberOfParametersForGene(0));
  }

  double getAicForGene(size_t i) const
  {
    return 2. * getNumberOfParametersForGene(i) - 2. * getLogLikelihoodForGene(i);
  }

  double getMeanForGene(size_t i, size_t partitionIndex) const
  {
    return distParamValues_[i][partitionIndex * 2];
  }

  double getNoiseForGene(size_t i, size_t partitionIndex) const
  {
    // This is 1/ shape
    return 1. / distParamValues_[i][partitionIndex * 2 + 1];
  }

  void printGeneParameters(ostream& out) const
  {
    out << "Gene,r";
    for (auto p : distParamNames_)
    {
      out << "," << p;
    }
    out << ",logLk,nbParam,AIC" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i] << "," << rs_[i];
      for (auto p : distParamValues_[i])
      {
        out << "," << p;
      }
      out << "," << getLogLikelihoodForGene(i) << "," << getNumberOfParametersForGene(i) << "," << getAicForGene(i) << endl;
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
    out << "Gene";
    for (auto& p : partition_.partitionMap)
    {
      out << "," << p.first << ".Mean," << p.first << ".Noise";
    }
    out << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i];
      for (size_t p = 0; p < partition_.size(); ++p)
      {
        out << "," << getMeanForGene(i, p) << "," << getNoiseForGene(i, p);
      }
      out << endl;
    }
  }


public:
  class GeneLikelihood :
    public virtual FunctionInterface,
    public AbstractParameterAliasable
  {
private:
    shared_ptr<CellPartitionFreeParameterSgeModel> sgeModel_;
    size_t geneId_;
    double r_;
    d1_array_dbl distParams_;

public:
    GeneLikelihood(shared_ptr<CellPartitionFreeParameterSgeModel> sgeModel, size_t geneId) :
      AbstractParameterAliasable(""),
      sgeModel_(sgeModel),
      geneId_(geneId),
      r_(0),
      distParams_(sgeModel->nbDistParams_)
    {
      addParameter_(new Parameter("r", sgeModel->rs_[geneId], Parameter::PROP_CONSTRAINT_EX));
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
      r_ = parameter("r").getValue();
      for (size_t j = 0; j < sgeModel_->nbDistParams_; ++j)
      {
        distParams_[j] = parameter(sgeModel_->distParamNames_[j]).getValue();
      }
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihoodForGene(geneId_, r_, distParams_);
    }
  };

  class CellLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<CellPartitionFreeParameterSgeModel> sgeModel_;
    size_t cellId_;
    double l_;

public:
    CellLikelihood(shared_ptr<CellPartitionFreeParameterSgeModel> sgeModel, size_t cellId) :
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

  void optimize(
      shared_ptr<OutputStream> profiler,
      shared_ptr<OutputStream> messenger,
      const std::vector<OptimizationParameters>& opt)
  {
    for (const auto& ppt : opt)
    {
      ApplicationTools::displayMessage("* " + ppt.description);

      // Cell-specific parameters:
      if (ppt.estimateNorms)
      {
        ApplicationTools::displayTask("  - Cell-specific parameters"); cout << endl;
        unsigned int progress = 0;
        for_each (std::execution::par,
                  cellIndex_.begin(),
                  cellIndex_.end(),
                  [&](size_t j)
        {
          auto f = make_shared<CellLikelihood>(shared_from_this(), j);
          auto optimizer = make_unique<BrentOneDimension>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->setInitialInterval(0.001, 100);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(ppt.tolerance);
          optimizer->optimize();
          ls_[j] = optimizer->getParameterValue("l");
          progress++;
          ApplicationTools::displayGauge(progress, cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("    Log likelihood is", getLogLikelihood());
      }
      if (ppt.outputCells != "none" && ppt.outputCells != "None")
      {
        ofstream streamCells(ppt.outputCells, ios::out);
        printCellParameters(streamCells);
        streamCells.close();
      }

      // Gene-specific parameters:
      ApplicationTools::displayTask("  - Gene distributions"); cout << endl;
      unsigned int progress = 0;
      for_each (std::execution::par,
                geneIndex_.begin(),
                geneIndex_.end(),
                [&](size_t i)
      {
        auto f = make_shared<GeneLikelihood>(shared_from_this(), i);
        if (ppt.aliasShape) {
          shapeIsAliased_ = true;
          auto pnames = f->getParameters().getParameterNames();
          pnames.erase(remove_if(pnames.begin(),
                                 pnames.end(),
                                 [](const string& name) {
                                   return(!name.ends_with("Shape"));
                                 }),
                       pnames.end());
          for (auto name : pnames) {
            if (name != pnames[0] ) {
              f->aliasParameters(pnames[0], name);
            }
          }
        } else {
          shapeIsAliased_ = true;
        }
        auto optimizer = make_unique<PowellMultiDimensions>(f);
        optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
        optimizer->setMessageHandler(messenger);
        optimizer->setProfiler(profiler);
        optimizer->setVerbose(0);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        optimizer->init(f->getIndependentParameters());
        optimizer->getStopCondition()->setTolerance(ppt.tolerance);
        optimizer->optimize();
        rs_[i] = optimizer->getParameterValue("r");
        for (size_t k = 0; k < nbDistParams_; ++k)
        {
          distParamValues_[i][k] = f->getParameterValue(distParamNames_[k]);
        }
        progress++;
        ApplicationTools::displayGauge(progress, geneIndex_.size() - 1, '=');
      });
      ApplicationTools::displayTaskDone();
      ApplicationTools::displayResult("    Log likelihood is", getLogLikelihood());
 
      if (ppt.outputGenes != "none" && ppt.outputGenes != "None")
      {
        ofstream streamGenes(ppt.outputGenes, ios::out);
        printGeneParameters(streamGenes);
        streamGenes.close();
      }
    }
  }

};

#endif //_CELL_PARTITION_FREE_PARAMETER_SGE_MODEL_H_

