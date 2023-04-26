/*
 * =====================================================================================
 *
 *       Filename:  nisce.cpp
 *
 *    Description:  Noise Inference from Single Cell Expression
 *
 *        Version:  1.0
 *        Created:  20.04.2023 10:39:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Julien Dutheil dutheil@evolbio.mpg.de
 *        Company:  Max Planck Institute for Evolutionary Biology
 *
 * =====================================================================================
 */


// Compile with:
// CC=gcc-12 CXX=g++-12 cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local .
// make
// Note: Clang not yet suported as C++17 does not implement parallelization.

// From the STL:
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <numeric>
#include <memory>
#include <cfloat>
#include <cmath>
#include <cctype>
#include <limits>
#include <execution>

// From Boost:
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

// From bpp-core:
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Exceptions.h>

using namespace std;
using namespace bpp;

using d2_array_ui = vector<vector<unsigned int> >;
using d1_array_ui = vector<unsigned int>;
using d3_array_dbl = vector<vector<vector<double> > >;
using d2_array_dbl = vector<vector<double> >;
using d1_array_dbl = vector<double>;
using d1_array_str = vector<string>;

// LogSumExp function, https://stackoverflow.com/questions/45943750/calculating-log-sum-exp-function-in-c
template<typename Iter>
typename iterator_traits<Iter>::value_type logSumExp(Iter begin, Iter end)
{
  using VT = typename iterator_traits<Iter>::value_type;
  if (begin == end)
    return VT{}
  ;
  auto max_elem = *max_element(begin, end);
  auto sum = accumulate(begin, end, VT{},
                        [max_elem](VT a, VT b) {
    return a + exp(b - max_elem);
  });
  return max_elem + log(sum);
}

double dpoisson(unsigned int x, double lambda)
{
  if (lambda == 0)
    return x == 0 ? 1.0 : 0.0;
  if (x < 0)
    return 0.;
  if (x <= lambda * DBL_MIN)
    return exp(-lambda);
  double p = exp(-lambda + x * log(lambda) - lgamma(x + 1));
  return max(p, std::numeric_limits< double >::min());
}

class SGEModel
{
public:
  d2_array_ui counts_;
  d1_array_str geneIds_;
  d1_array_str cellIds_;
  size_t nbGenes_;
  size_t nbCells_;
  string model_;
  size_t nbCategories_;

  // Distribution parameters
  size_t nbDistParams_;
  vector<string> distParamNames_;
  vector<shared_ptr<Constraint> > distParamConstraints_;
  d2_array_dbl distParamValues_;

  d1_array_dbl rs_;           // Dropout parameters
  d1_array_dbl ls_;           // Library size parameters
  d3_array_dbl distArray_;    // Intermediate calculations. Contains, for each gene, each category, two numbers: value and probability.

public:
  SGEModel(
    const d2_array_ui& counts,
    const d1_array_str& geneIds,
    const d1_array_str& cellIds,
    const string& model,
    size_t nbCategories) :
    counts_(counts),
    geneIds_(geneIds),
    cellIds_(cellIds),
    nbGenes_(counts.size()),
    nbCells_(nbGenes_ == 0 ? 0 : counts[0].size()),
    model_(model),
    nbCategories_(nbCategories),
    nbDistParams_(0),  // Will be initialized in the constructor body
    distParamNames_(), // Will be initialized in the constructor body
    distParamConstraints_(), // Will be initialized in the constructor body
    distParamValues_(nbGenes_),
    rs_(nbGenes_),
    ls_(nbCells_),
    distArray_(nbGenes_)
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
    // Dropout probability is proportion of 0s
    // First distribution parameter is the mean
    for (size_t i  = 0; i < nbGenes_; ++i)
    {
      double c = static_cast<double>(count(counts_[i].begin(), counts_[i].end(), 0));
      double s = accumulate(counts_[i].begin(), counts_[i].end(), 0.);
      rs_[i] = c / static_cast<double>(nbCells_);
      distParamValues_[i][0] = s / (static_cast<double>(nbCells_) - c);
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
  }

  /**
   * @brief Output information about the model
   */
  void print(ostream& out) const
  {
    out << "-----------------------------------------------------------------" << endl;
    out << " SGE model with " << nbGenes_ << " genes and " << nbCells_ << " cells." << endl;
    out << " Poisson distribution is modulated by a " << model_ << " distribution." << endl;
    out << "-----------------------------------------------------------------" << endl;
  }

  /**
   * @brief Likelihood calculations
   */
  double computeLogLikelihoodForGene(size_t i, double r, const d1_array_dbl& distParams)
  {
    computeCategories(i, distParams);
    double logLikelihood = 0.;
    for (size_t j = 0; j < nbCells_; ++j)
    {
      d1_array_dbl s(nbCategories_);
      for (size_t k = 0; k < nbCategories_; ++k)
      {
        s[k] = dpoisson(counts_[i][j], distArray_[i][k][0] / ls_[j]) * distArray_[i][k][1];
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
    if (isnan(logLikelihood) || isinf(logLikelihood))
    {
      logLikelihood = 0.;
      cout << "Gene " << i << endl;
      for (size_t j = 0; j < nbCells_; ++j)
      {
        cout << "Cell " << j << endl;
        d1_array_dbl s(nbCategories_);
        for (size_t k = 0; k < nbCategories_; ++k)
        {
          s[k] = dpoisson(counts_[i][j], distArray_[i][k][0] / ls_[j]) * distArray_[i][k][1];
          cout << k << ", counts = " << counts_[i][j] << ", distArray = " << distArray_[i][k][0] <<  ", dpoisson = " << dpoisson(counts_[i][j], distArray_[i][k][0] / ls_[j]) << endl;
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
        cout << i << "," << j << ", S = " << S << ", P = " << P << ", log(P) = " << log(P) << endl;
        logLikelihood += log(P);
      }
      cout << logLikelihood << endl;

      throw (exception());
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
        s[k] = dpoisson(counts_[i][j], distArray_[i][k][0] / l) * distArray_[i][k][1];
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
      logLikelihood = 0.;
      cout << "Cell " << j << endl;
      for (size_t i = 0; i < nbGenes_; ++i)
      {
        cout << "Gene " << i << endl;
        d1_array_dbl s(nbCategories_);
        for (size_t k = 0; k < nbCategories_; ++k)
        {
          s[k] = dpoisson(counts_[i][j], distArray_[i][k][0] / l) * distArray_[i][k][1];
          cout << k << ", counts = " << counts_[i][j] << ", disrArray = " << distArray_[i][k][0] <<  ", dpoisson = " << dpoisson(counts_[i][j], distArray_[i][k][0] / l) << endl;
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
        cout << i << "," << j << ", S = " << S << ", P = " << P << ", log(P) = " << log(P) << endl;
        logLikelihood += log(P);
      }
      cout << logLikelihood << endl;

      throw (exception());
    }
    return logLikelihood;
  }

  double getLogLikelihoodForGene(size_t i) const {
    return const_cast<SGEModel*>(this)->computeLogLikelihoodForGene(i, rs_[i], distParamValues_[i]); //Actually ok to do so here, since we are not modifying the object.
  }

  double getAicForGene(size_t i) const {
    return 2. * static_cast<double>(nbDistParams_ + 1) - 2. * getLogLikelihoodForGene(i); //Needs to account for r too (hence the +1).
  }

  double getMeanForGene(size_t i) const {
    if (model_ == "Gamma") {
      //This is the mean parameter:
      return distParamValues_[i][0];
    } else if (model_ == "ScaledBeta") {
      //This is scale * alpha / (alpha + beta)
      return distParamValues_[i][0] * distParamNames_[i][1] / (distParamValues_[i][1] + distParamValues_[i][2]);
    }
    throw Exception("Unsupported model: " + model_);
  }

  double getNoiseForGene(size_t i) const {
    if (model_ == "Gamma") {
      //This is 1/ shape
      return 1. / distParamValues_[i][1];
    } else if (model_ == "ScaledBeta") {
      //This is beta / (alpha * (alpha + beta + 1))
      return distParamNames_[i][2] / (distParamValues_[i][1] * (distParamValues_[i][1] + distParamValues_[i][2] + 1.));
    }
    throw Exception("Unsupported model: " + model_);
  }

  void printGeneParameters(ostream& out) const
  {
    out << "Gene,r";
    for (auto p : distParamNames_)
      out << "," << p;
    out << ",logLk,nbParam,AIC" << endl;

    for (size_t i = 0; i < nbGenes_; ++i)
    {
      out << geneIds_[i] << "," << rs_[i];
      for (auto p : distParamValues_[i])
        out << "," << p;
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
};


class GeneLikelihood :
  public virtual FunctionInterface,
  public AbstractParametrizable
{
private:
  shared_ptr<SGEModel> sgeModel_;
  size_t geneId_;
  double r_;
  d1_array_dbl distParams_;

public:
  GeneLikelihood(shared_ptr<SGEModel> sgeModel, size_t geneId) :
    AbstractParametrizable(""),
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
    r_ = getParameter("r").getValue();
    for (size_t j = 0; j < sgeModel_->nbDistParams_; ++j)
    {
      distParams_[j] = getParameter(sgeModel_->distParamNames_[j]).getValue();
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
  shared_ptr<SGEModel> sgeModel_;
  size_t cellId_;
  double l_;

public:
  CellLikelihood(shared_ptr<SGEModel> sgeModel, size_t cellId) :
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

d2_array_ui readCounts(istream& input, d1_array_str& geneIds, d1_array_str& cellIds, char sep = ',')
{
  geneIds.resize(0);
  cellIds.resize(0);
  d2_array_ui counts;
  string line;
  getline(input, line); // Skip header
  stringstream headerStream(line);
  string value;
  getline(headerStream, value, sep); //Remove gene id col header (might just be empty).
  while (getline(headerStream, value, sep))
  {
    cellIds.push_back(value);
  }
  size_t nbCells = cellIds.size();
  while (getline(input, line))
  {
    stringstream lineStream(line);
    d1_array_ui row;
    getline(lineStream, value, sep);
    geneIds.push_back(value); // Add gene id
    while (getline(lineStream, value, sep))
    {
      stringstream valueStream(value);
      unsigned int count;
      valueStream >> count;
      row.push_back(count);
    }
    if (row.size() == nbCells)
    {
      counts.push_back(row);
    }
    else
    {
      cerr << "Error, input file has unequal number of rows." << endl;
      throw exception();
    }
    ApplicationTools::displayUnlimitedGauge(counts.size(), "Parsing input file: ");
  }
  return counts;
}

void help()
{
  cout << "This is NISCE 0.01." << endl;
  cout << "Usage:" << endl;
  cout << "  nisce expression_file.csv" << endl;
}


int main(int argc, char* argv[])
{
  cout << "******************************************************************" << endl;
  cout << "*                  N.I.S.C.E., version 0.0.1                     *" << endl;
  cout << "* Author: J. Dutheil                        Created on  24/04/23 *" << endl;
  cout << "*                                           Last Modif. 24/04/23 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (argc == 1)
  {
    help();
    exit(0);
  }

  try
  {
    BppApplication nisce(argc, argv, "nisce");
    nisce.startTimer();

    string inputPath = ApplicationTools::getAFilePath("input.file", nisce.getParams(), true, true);
    char sep = ApplicationTools::getParameter<char>("input.file.column_separator", nisce.getParams(), ',');

    // Read the counts:
    ifstream inputFile(inputPath, ios::in);
    vector<string> geneIds;
    vector<string> cellIds;
    auto counts = readCounts(inputFile, geneIds, cellIds, sep);
    if (counts.size() == 0)
    {
      cout << "Error, file " << argv[1] << " seems to contain no data." << endl;
      throw exception();
    }
    ApplicationTools::displayTaskDone();
    ApplicationTools::displayResult("Number of genes in file", counts.size());
    ApplicationTools::displayResult("Number of cells in file", counts[0].size());

    // Discard empty genes:
    ApplicationTools::displayTask("Filter 'empty' genes");
    size_t minGeneObs = ApplicationTools::getParameter<size_t>("min_nb_cells_with_obs", nisce.getParams(), 1);
    if (minGeneObs < 1)
      throw Exception("Minimum number of observations per gene must be at least 1.");
    unsigned int geneCount = 0;
    for (size_t i = counts.size(); i > 0; --i)
    {
      size_t s = counts[0].size() - count(counts[i - 1].begin(), counts[i - 1].end(), 0);
      if (s < minGeneObs)
      {
        // cout << "Gene " << i << " has " << s << " counts only and was discarded." << endl;
        counts.erase(counts.begin() + static_cast<ptrdiff_t>(i) - 1, counts.begin() + static_cast<ptrdiff_t>(i));
        geneIds.erase(geneIds.begin() + static_cast<ptrdiff_t>(i) - 1, geneIds.begin() + static_cast<ptrdiff_t>(i));
        geneCount++;
      }
    }
    ApplicationTools::displayTaskDone();
    ApplicationTools::displayResult("Number of discarded genes", geneCount);

    // Discard empty cells:
    ApplicationTools::displayTask("Filter 'empty' cells");
    vector<ptrdiff_t> cellsToDiscard;
    for (size_t j = 0; j < counts[0].size(); ++j)
    {
      size_t i;
      for (i = 0; i < counts.size() && counts[i][j] == 0; ++i) {}
      if (i == counts.size())
        cellsToDiscard.push_back(static_cast<ptrdiff_t>(j));
    }
    for (auto& v : counts)
    {
      for (auto it = cellsToDiscard.rbegin(); it != cellsToDiscard.rend(); ++it)
      {
        v.erase(v.begin() + *it, v.begin() + *it + 1);
        cellIds.erase(cellIds.begin() + *it, cellIds.begin() + *it + 1);
      }
    }
    ApplicationTools::displayTaskDone();
    ApplicationTools::displayResult("Number of discarded cells", cellsToDiscard.size());

    // Initialize parameters:
    string model = ApplicationTools::getStringParameter("gene_distribution", nisce.getParams(), "ScaledBeta");
    size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 5);

    // Create the data structure:
    auto sgeModel = make_shared<SGEModel>(counts, geneIds, cellIds, model, nbCategories);
    sgeModel->print(cout);
    sgeModel->initParameters();

    // Now optimize parameters.

    double tol1 = ApplicationTools::getDoubleParameter("tolerance1", nisce.getParams(), 0.01);

    string prPath = ApplicationTools::getAFilePath("profiler", nisce.getParams(), false, false);
    string mhPath = ApplicationTools::getAFilePath("message_handler", nisce.getParams(), false, false);

    auto profiler =
      (prPath == "none") ? 0 :
      (prPath == "std") ? ApplicationTools::message :
      make_shared<StlOutputStream>(make_unique<fstream>(prPath.c_str(), ios::out));
    if (profiler)
      profiler->setPrecision(20);
    ApplicationTools::displayResult("-- Optimization profile in", prPath);

    auto messenger =
      (mhPath == "none") ? nullptr :
      (mhPath == "std") ? ApplicationTools::message :
      make_shared<StlOutputStream>(make_unique<ofstream>(mhPath.c_str(), ios::out));
    if (messenger)
      messenger->setPrecision(20);
    ApplicationTools::displayResult("-- Optimization messages in", mhPath);

    // First gene specific parameters:
    vector<size_t> genes;
    for (size_t i = 0; i < sgeModel->nbGenes_; ++i)
    {
      genes.push_back(i);
    }

    ApplicationTools::displayTask("Coarse fitting of gene distributions");
    cout << endl;
    unsigned int progress = 0;
    for_each (std::execution::par, genes.begin(), genes.end(), [&](size_t i) {
      auto f = make_shared<GeneLikelihood>(sgeModel, i);
      auto optimizer = make_unique<PowellMultiDimensions>(f);
      optimizer->setMessageHandler(messenger);
      optimizer->setProfiler(profiler);
      optimizer->setVerbose(0);
      optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      optimizer->init(f->getParameters());
      optimizer->getStopCondition()->setTolerance(tol1);
      optimizer->optimize();
      sgeModel->rs_[i] = optimizer->getParameterValue("r");
      for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
      {
        sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
      }
      progress++;
      ApplicationTools::displayGauge(progress, genes.size() - 1, '=');
    });
    ApplicationTools::displayTaskDone();

    string outputGenes1 = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
    if (outputGenes1 != "none") {
      ofstream streamGenes1(outputGenes1, ios::out);
      sgeModel->printGeneParameters(streamGenes1);
      streamGenes1.close();
    }

    // Then cell specific parameters:
    vector<size_t> cells;
    for (size_t j = 0; j < sgeModel->nbCells_; ++j)
    {
      cells.push_back(j);
    }

    double tol2 = ApplicationTools::getDoubleParameter("tolerance2", nisce.getParams(), 0.0001);

    ApplicationTools::displayTask("Fitting of cell-specific parameters");
    cout << endl;
    progress = 0;
    for_each (std::execution::par, cells.begin(), cells.end(), [&](size_t j) {
      auto f = make_shared<CellLikelihood>(sgeModel, j);
      auto optimizer = make_unique<BrentOneDimension>(f);
      optimizer->setMessageHandler(messenger);
      optimizer->setProfiler(profiler);
      optimizer->setVerbose(0);
      optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      optimizer->setInitialInterval(0.001, 100);
      optimizer->init(f->getParameters());
      optimizer->getStopCondition()->setTolerance(tol2);
      optimizer->optimize();
      sgeModel->ls_[j] = optimizer->getParameterValue("l");
      progress++;
      ApplicationTools::displayGauge(progress, cells.size() - 1, '=');
    });
    ApplicationTools::displayTaskDone();

    string outputCells1 = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
    if (outputCells1 != "none") {
      ofstream streamCells1(outputCells1, ios::out);
      sgeModel->printCellParameters(streamCells1);
      streamCells1.close();
    }

    // Re-fit distribution parameters
    double tol3 = ApplicationTools::getDoubleParameter("tolerance3", nisce.getParams(), 0.000001);
    
    ApplicationTools::displayTask("Precise fitting of gene distributions");
    cout << endl;
    progress = 0;
    for_each (std::execution::par, genes.begin(), genes.end(), [&](size_t i) {
      auto f = make_shared<GeneLikelihood>(sgeModel, i);
      auto optimizer = make_unique<PowellMultiDimensions>(f);
      optimizer->setMessageHandler(messenger);
      optimizer->setProfiler(profiler);
      optimizer->setVerbose(0);
      optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      optimizer->init(f->getParameters());
      optimizer->getStopCondition()->setTolerance(tol3);
      optimizer->optimize();
      sgeModel->rs_[i] = optimizer->getParameterValue("r");
      for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
      {
        sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
      }
      progress++;
      ApplicationTools::displayGauge(progress, genes.size() - 1, '=');
    });
    ApplicationTools::displayTaskDone();

    string outputGenes2 = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
    if (outputGenes2 != "none") {
      ofstream streamGenes2(outputGenes2, ios::out);
      sgeModel->printGeneParameters(streamGenes2);
      streamGenes2.close();
    }

    string outputGeneResult = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
    if (outputGenes2 != "none") {
      ofstream streamGeneResults(outputGeneResult, ios::out);
      sgeModel->printResultsForGenes(streamGeneResults);
      streamGeneResults.close();
    }

    nisce.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return 0;
}

//TODO
//- in final step, compute both gamma and beta and keep for each gene the one with lowest AIC?? => output logLik per gene!
