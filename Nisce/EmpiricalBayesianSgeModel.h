/*
 * =====================================================================================
 *
 *       Filename:  EmpiricalBayesianSgeModel.h
 *
 * =====================================================================================
 */

#include "ProbFunctions.h"

// Only supports Gamma model. We assume all counts to share independent Gammas, but we could inforce a single one (like ASRV vs SSRV phylogenetic models).
class EmpiricalBayesianSgeModel
{
public:
  d2_array_ui counts_;
  d1_array_str geneIds_;
  d1_array_str cellIds_;
  size_t nbGenes_;
  size_t nbCells_;
  size_t nbMeanCategories_;
  size_t nbShapeCategories_;
  size_t nbCategories_;

  // Distribution parameters
  double meanMean_;
  double shapeMean_;
  double meanShape_;
  double shapeShape_;

  d1_array_dbl distMeanValues_;
  d1_array_dbl distShapeValues_;

  d1_array_dbl rs_;           // Dropout parameters
  d1_array_dbl ls_;           // Library size parameters
  d3_array_dbl distArray_;    // Intermediate calculations. Contains, for each param category, each lambda category, two numbers: value and probability.

  vector<size_t> geneIndex_;

public:
  EmpiricalBayesianSgeModel(
    const d2_array_ui& counts,
    const d1_array_str& geneIds,
    const d1_array_str& cellIds,
    size_t nbMeanCategories,
    size_t nbShapeCategories,
    size_t nbCategories) :
    counts_(counts),
    geneIds_(geneIds),
    cellIds_(cellIds),
    nbGenes_(counts.size()),
    nbCells_(nbGenes_ == 0 ? 0 : counts[0].size()),
    nbMeanCategories_(nbMeanCategories),
    nbShapeCategories_(nbShapeCategories),
    nbCategories_(nbCategories),
    meanMean_(1.),
    shapeMean_(1.),
    meanShape_(1.),
    shapeShape_(1.),
    distMeanValues_(nbMeanCategories),
    distShapeValues_(nbShapeCategories),
    rs_(nbGenes_),
    ls_(nbCells_),
    distArray_(nbMeanCategories * nbShapeCategories),
    geneIndex_()
  {
    for (auto& a : distArray_)
    {
      a.resize(nbCategories_);
      for (auto& b : a)
      {
        b.resize(2);
      }
    }
    
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      geneIndex_.push_back(i);
    }
  }

public:
  /**
   * @brief Recomputes the Discretized Poisson modulated PDF for a given gene.
   */
  void computeCategories()
  {
    for (size_t i = 0; i < nbMeanCategories_; ++i)
    {
      distMeanValues_[i] = meanMean_ * (1. / shapeMean_) * boost::math::gamma_p_inv(shapeMean_, ((double)i + 0.5) / (double)nbCategories_); // We use the median as category value
    }
    for (size_t i = 0; i < nbShapeCategories_; ++i)
    {
      distShapeValues_[i] = meanShape_ * (1. / shapeShape_) * boost::math::gamma_p_inv(shapeShape_, ((double)i + 0.5) / (double)nbCategories_); // We use the median as category value
    }
    for (size_t k = 0; k < nbMeanCategories_; ++k)
    {
      for (size_t j = 0; j < nbShapeCategories_; ++j)
      {
        for (size_t i = 0; i < nbCategories_; ++i)
        {
          distArray_[k * nbMeanCategories_ + j][i][0] = distMeanValues_[k] * (1. / distShapeValues_[j]) * boost::math::gamma_p_inv(distShapeValues_[j], ((double)i + 0.5) / (double)nbCategories_); // We use the median as category value
          distArray_[k * nbMeanCategories_ + j][i][1] = 1. / static_cast<double>(nbCategories_);
        }
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
      double c = static_cast<double>(count(counts_[i].begin(), counts_[i].end(), 0));
      rs_[i] = c / static_cast<double>(nbCells_);
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
    out << " Gamma param. are estimated via an empirical Bayesian approach." << endl;
    out << "-----------------------------------------------------------------" << endl;
  }

  /**
   * @brief Likelihood calculations
   */
  double computeLogLikelihoodForGene(size_t i, double r)
  {
    double logLikelihood = 0.;
    for (size_t j = 0; j < nbCells_; ++j)
    {
      double l = ls_[j];
      d1_array_dbl s(nbMeanCategories_ * nbShapeCategories_ * nbCategories_);
      size_t z = 0;
      for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
      {
        for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
        {
          for (size_t k3 = 0; k3 < nbCategories_; ++k3)
          {
            s[z++] = dpoisson(counts_[i][j], distArray_[k1 * nbMeanCategories_ + k2][k3][0] * l) * distArray_[k1 * nbMeanCategories_ + k2][k3][1];
          }
        }
      }
      auto S = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbMeanCategories_ * nbShapeCategories_);
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
      d1_array_dbl s(nbMeanCategories_ * nbShapeCategories_ * nbCategories_);
      size_t z = 0;
      for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
      {
        for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
        {
          for (size_t k3 = 0; k3 < nbCategories_; ++k3)
          {
            s[z++] = dpoisson(counts_[i][j], distArray_[k1 * nbMeanCategories_ + k2][k3][0] * l) * distArray_[k1 * nbMeanCategories_ + k2][k3][1];
          }
        }
      }
      auto S = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbMeanCategories_ * nbShapeCategories_);
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

  double computeLogLikelihood(double meanMean, double shapeMean, double meanShape, double shapeShape)
  {
    meanMean_ = meanMean;
    shapeMean_ = shapeMean;
    meanShape_ = meanShape;
    shapeShape_ = shapeShape;
    computeCategories();
    d1_array_dbl logLik(nbGenes_);

    for_each (std::execution::par, geneIndex_.begin(), geneIndex_.end(), [&](size_t i) {
      logLik[i] = computeLogLikelihoodForGene(i, rs_[i]);
    });
    sort(logLik.begin(), logLik.end()); // Sort by increasing value to limit underflow.
    double logLikelihood = accumulate(logLik.begin(), logLik.end(), 0.);

    return logLikelihood;
  }

  double getLogLikelihoodForGene(size_t i) const
  {
    return const_cast<EmpiricalBayesianSgeModel*>(this)->computeLogLikelihoodForGene(i, rs_[i]); // Actually ok to do so here, since we are not modifying the object.
  }

  double getMeanForGeneForCell(size_t i, size_t j) const
  {
    double l = ls_[j];
    d1_array_dbl condLik(nbMeanCategories_);
    for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
    {
      d1_array_dbl s(nbCategories_ * nbShapeCategories_);
      size_t z = 0;
      for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
      {
        for (size_t k3 = 0; k3 < nbCategories_; ++k3)
        {
          s[z++] = dpoisson(counts_[i][j], distArray_[k1 * nbMeanCategories_ + k2][k3][0] * l) * distArray_[k1 * nbMeanCategories_ + k2][k3][1] ;
        }
      }
      condLik[k1] = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbMeanCategories_ + nbShapeCategories_);
    }
    double probSum = accumulate(condLik.begin(), condLik.end(), 0.);

    double postMean = 0.;
    for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
    {
      postMean += distMeanValues_[k1] * condLik[k1] / probSum; 
    }
    return postMean;
  }

  double getMeanEntropyForGeneForCell(size_t i, size_t j) const
  {
    double l = ls_[j];
    d1_array_dbl condLik(nbMeanCategories_);
    for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
    {
      d1_array_dbl s(nbCategories_ * nbShapeCategories_);
      size_t z = 0;
      for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
      {
        for (size_t k3 = 0; k3 < nbCategories_; ++k3)
        {
          s[z++] = dpoisson(counts_[i][j], distArray_[k1 * nbMeanCategories_ + k2][k3][0] * l) * distArray_[k1 * nbMeanCategories_ + k2][k3][1] ;
        }
      }
      condLik[k1] = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbMeanCategories_ + nbShapeCategories_);
    }
    double probSum = accumulate(condLik.begin(), condLik.end(), 0.);

    double ent = 0.;
    for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
    {
      double p = condLik[k1] / probSum; 
      ent -= p * log(p); 
    }
    return ent;
  }

  double getNoiseForGeneForCell(size_t i, size_t j) const
  {
    double l = ls_[j];
    d1_array_dbl condLik(nbShapeCategories_);
    for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
    {
      d1_array_dbl s(nbCategories_ * nbMeanCategories_);
      size_t z = 0;
      for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
      {
        for (size_t k3 = 0; k3 < nbCategories_; ++k3)
        {
          s[z++] = dpoisson(counts_[i][j], distArray_[k1 * nbMeanCategories_ + k2][k3][0] * l) * distArray_[k1 * nbMeanCategories_ + k2][k3][1] ;
        }
      }
      condLik[k2] = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbMeanCategories_ + nbShapeCategories_);
    }
    double probSum = accumulate(condLik.begin(), condLik.end(), 0.);

    double postNoise = 0.;
    for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
    {
      postNoise += (1. / distShapeValues_[k2]) * condLik[k2] / probSum; 
    }
    return postNoise;
  }

  double getNoiseEntropyForGeneForCell(size_t i, size_t j) const
  {
    double l = ls_[j];
    d1_array_dbl condLik(nbShapeCategories_);
    for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
    {
      d1_array_dbl s(nbCategories_ * nbMeanCategories_);
      size_t z = 0;
      for (size_t k1 = 0; k1 < nbMeanCategories_; ++k1)
      {
        for (size_t k3 = 0; k3 < nbCategories_; ++k3)
        {
          s[z++] = dpoisson(counts_[i][j], distArray_[k1 * nbMeanCategories_ + k2][k3][0] * l) * distArray_[k1 * nbMeanCategories_ + k2][k3][1] ;
        }
      }
      condLik[k2] = accumulate(s.begin(), s.end(), 0.) / static_cast<double>(nbMeanCategories_ + nbShapeCategories_);
    }
    double probSum = accumulate(condLik.begin(), condLik.end(), 0.);

    double ent = 0.;
    for (size_t k2 = 0; k2 < nbShapeCategories_; ++k2)
    {
      double p = condLik[k2] / probSum; 
      ent -= p * log(p); 
    }
    return ent;
  }


  void printHyperParameters(ostream& out) const
  {
    out << "- Mean hyper parameters:" << endl;
    out << "  * Mean: " << meanMean_ << endl;
    out << "  * Shape: " << shapeMean_ << endl;
    out << "- Shape hyper parameters:" << endl;
    out << "  * Mean: " << meanShape_ << endl;
    out << "  * Shape: " << shapeShape_ << endl;
  }

  void printGeneParameters(ostream& out) const
  {
    out << "Gene,r,logLk" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i] << "," << rs_[i] << "," << getLogLikelihoodForGene(i) << endl;
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

  void printMeanResults(ostream& out) const
  {
    out << "Gene";
    for (size_t j = 0; j < nbCells_; ++j)
    {
      out << "," << cellIds_[j];
    }
    out << endl;
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i];
      for (size_t j = 0; j < nbCells_; ++j)
      {
	out << "," << getMeanForGeneForCell(i, j);
      }
      out << endl;
    }
  }

  void printNoiseResults(ostream& out) const
  {
    out << "Gene";
    for (size_t j = 0; j < nbCells_; ++j)
    {
      out << "," << cellIds_[j];
    }
    out << endl;
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i];
      for (size_t j = 0; j < nbCells_; ++j)
      {
	out << "," << getNoiseForGeneForCell(i, j);
      }
      out << endl;
    }
  }

  void printMeanEntropyResults(ostream& out) const
  {
    out << "Gene";
    for (size_t j = 0; j < nbCells_; ++j)
    {
      out << "," << cellIds_[j];
    }
    out << endl;
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i];
      for (size_t j = 0; j < nbCells_; ++j)
      {
	out << "," << getMeanEntropyForGeneForCell(i, j);
      }
      out << endl;
    }
  }

  void printNoiseEntropyResults(ostream& out) const
  {
    out << "Gene";
    for (size_t j = 0; j < nbCells_; ++j)
    {
      out << "," << cellIds_[j];
    }
    out << endl;
    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i];
      for (size_t j = 0; j < nbCells_; ++j)
      {
	out << "," << getNoiseEntropyForGeneForCell(i, j);
      }
      out << endl;
    }
  }

public:
  class EBLikelihood :
    public virtual FunctionInterface,
    public AbstractParametrizable
  {
private:
    shared_ptr<EmpiricalBayesianSgeModel> sgeModel_;
    double meanMean_;
    double shapeMean_;
    double meanShape_;
    double shapeShape_;

public:
    EBLikelihood(shared_ptr<EmpiricalBayesianSgeModel> sgeModel) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      meanMean_(0),
      shapeMean_(0),
      meanShape_(0),
      shapeShape_(0)
    {
      addParameter_(new Parameter("meanMean", sgeModel->meanMean_, make_shared<IntervalConstraint>(0.01, 999, true, true)));
      addParameter_(new Parameter("shapeMean", sgeModel->shapeMean_, make_shared<IntervalConstraint>(0.01, 999, true, true)));
      addParameter_(new Parameter("meanShape", sgeModel->meanShape_, make_shared<IntervalConstraint>(0.01, 999, true, true)));
      addParameter_(new Parameter("shapeShape", sgeModel->shapeShape_, make_shared<IntervalConstraint>(0.01, 999, true, true)));
    }

    EBLikelihood* clone() const override { return new EBLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      meanMean_ = getParameter("meanMean").getValue();
      shapeMean_ = getParameter("shapeMean").getValue();
      meanShape_ = getParameter("meanShape").getValue();
      shapeShape_ = getParameter("shapeShape").getValue();
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihood(meanMean_, shapeMean_, meanShape_, shapeShape_);
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

public:
    GeneLikelihood(shared_ptr<EmpiricalBayesianSgeModel> sgeModel, size_t geneId) :
      AbstractParametrizable(""),
      sgeModel_(sgeModel),
      geneId_(geneId),
      r_(0)
    {
      addParameter_(new Parameter("r", sgeModel->rs_[geneId], Parameter::PROP_CONSTRAINT_EX));
    }

    GeneLikelihood* clone() const override { return new GeneLikelihood(*this); }

public:
    void setParameters(const ParameterList& parameters) override
    {
      matchParametersValues(parameters);
      r_ = getParameter("r").getValue();
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihoodForGene(geneId_, r_);
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
      l_ = getParameter("l").getValue();
    }

    double getValue() const override
    {
      return -sgeModel_->computeLogLikelihoodForCell(cellId_, l_);
    }
  };
};


