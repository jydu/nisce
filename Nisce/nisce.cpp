// SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * =====================================================================================
 *
 *       Filename:  nisce.cpp
 *
 *    Description:  Noise Inference from Single Cell Expression
 *
 *        Version:  0.0.1
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
// CC=gcc-12 CXX=g++-12 cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
// make
// Note: Clang not yet suported as C++17 does not implement parallelization.

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1 // _Suppress TBB warning

// From the STL:
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <memory>
#include <cfloat>
#include <cmath>
#include <cctype>
#include <limits>
#include <execution>

// From bpp-core:
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Exceptions.h>

using namespace std;
using namespace bpp;

#include "ProbFunctions.h"
#include "FreeParameterSgeModel.h"
#include "LessFreeParameterSgeModel.h"
#include "EmpiricalBayesianSgeModel.h"
#include "FullEmpiricalBayesianSgeModel.h"
#include "CellPartitionFreeParameterSgeModel.h"


d2_array_ui readCounts(istream& input, d1_array_str& geneIds, d1_array_str& cellIds, char sep = ',')
{
  geneIds.resize(0);
  cellIds.resize(0);
  d2_array_ui counts;
  string line;
  getline(input, line); // Skip header
  stringstream headerStream(line);
  string value;
  getline(headerStream, value, sep); // Remove gene id col header (might just be empty).
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
  cout << "This is NISCE 0.0.2." << endl;
  cout << "Usage:" << endl;
  cout << "  nisce param=options.bpp" << endl;
}


int main(int argc, char* argv[])
{
  cout << "******************************************************************" << endl;
  cout << "*                  N.I.S.C.E., version 0.0.3                     *" << endl;
  cout << "* Author: J. Dutheil                        Created on  24/04/23 *" << endl;
  cout << "*                                           Last Modif. 29/08/23 *" << endl;
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
      for (i = 0; i < counts.size() && counts[i][j] == 0; ++i)
      {}
      if (i == counts.size())
        cellsToDiscard.push_back(static_cast<ptrdiff_t>(j));
    }
    for (auto& v : counts)
    {
      for (auto it = cellsToDiscard.rbegin(); it != cellsToDiscard.rend(); ++it)
      {
        v.erase(v.begin() + *it, v.begin() + *it + 1);
      }
    }
    for (auto it = cellsToDiscard.rbegin(); it != cellsToDiscard.rend(); ++it)
    {
      cellIds.erase(cellIds.begin() + *it, cellIds.begin() + *it + 1);
    }
    ApplicationTools::displayTaskDone();
    ApplicationTools::displayResult("Number of discarded cells", cellsToDiscard.size());

    // Select method:
    string partitionFile = ApplicationTools::getStringParameter("cell_partition.file", nisce.getParams(), "None");
    string method = ApplicationTools::getStringParameter("method", nisce.getParams(), "FP");

    string prPath = ApplicationTools::getAFilePath("profiler", nisce.getParams(), false, false);
    string mhPath = ApplicationTools::getAFilePath("message_handler", nisce.getParams(), false, false);

    auto profiler =
      (prPath == "none") ? nullptr :
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

    if (partitionFile == "None" || partitionFile == "none")
    {
      /************************************
      *                                  *
      * Fit a single model to all cells. *
      *                                  *
      ************************************/

      if (method == "FP" || method == "FreeParameters")
      {
        /***********************************************************
        *                                                         *
        * Simple model with free parameters. Gamma or ScaledBeta. *
        *                                                         *
        ***********************************************************/

        // Initialize parameters:
        string model = ApplicationTools::getStringParameter("gene_distribution", nisce.getParams(), "ScaledBeta");
        size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 5);

        // Create the data structure:
        auto sgeModel = make_shared<FreeParameterSgeModel>(counts, geneIds, cellIds, model, nbCategories);
        sgeModel->print(cout);
        sgeModel->initParameters();

        // Now optimize parameters.
        double tol1 = ApplicationTools::getDoubleParameter("tolerance1", nisce.getParams(), 0.01);

        ApplicationTools::displayTask("Coarse fitting of gene distributions");
        cout << endl;
        unsigned int progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<FreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
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
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes1 = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
        if (outputGenes1 != "none")
        {
          ofstream streamGenes1(outputGenes1, ios::out);
          sgeModel->printGeneParameters(streamGenes1);
          streamGenes1.close();
        }

        // Then cell specific parameters:

        double tol2 = ApplicationTools::getDoubleParameter("tolerance2", nisce.getParams(), 0.0001);

        ApplicationTools::displayTask("Fitting of cell-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->cellIndex_.begin(),
                  sgeModel->cellIndex_.end(),
                  [&](size_t j)
        {
          auto f = make_shared<FreeParameterSgeModel::CellLikelihood>(sgeModel, j);
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
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputCells2 = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
        if (outputCells2 != "none")
        {
          ofstream streamCells2(outputCells2, ios::out);
          sgeModel->printCellParameters(streamCells2);
          streamCells2.close();
        }

        ApplicationTools::displayTask("Re-fitting of gene distributions");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<FreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->rs_[i] = optimizer->getParameterValue("r");
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes2 = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
        if (outputGenes2 != "none")
        {
          ofstream streamGenes2(outputGenes2, ios::out);
          sgeModel->printGeneParameters(streamGenes2);
          streamGenes2.close();
        }

        // Fine fitting
        
        double tol3 = ApplicationTools::getDoubleParameter("tolerance3", nisce.getParams(), 0.000001);

        ApplicationTools::displayTask("Fine-fitting of cell-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->cellIndex_.begin(),
                  sgeModel->cellIndex_.end(),
                  [&](size_t j) 
        {
          auto f = make_shared<FreeParameterSgeModel::CellLikelihood>(sgeModel, j);
          auto optimizer = make_unique<BrentOneDimension>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->setInitialInterval(0.001, 100);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol3);
          optimizer->optimize();
          sgeModel->ls_[j] = optimizer->getParameterValue("l");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputCells3 = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
        if (outputCells3 != "none")
        {
          ofstream streamCells3(outputCells3, ios::out);
          sgeModel->printCellParameters(streamCells3);
          streamCells3.close();
        }

        ApplicationTools::displayTask("Fine-fitting of gene distributions");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<FreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->rs_[i] = optimizer->getParameterValue("r");
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes3 = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
        if (outputGenes3 != "none")
        {
          ofstream streamGenes3(outputGenes3, ios::out);
          sgeModel->printGeneParameters(streamGenes3);
          streamGenes3.close();
        }

        // Final results:
        string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
        if (outputGeneResults != "none")
        {
          ofstream streamGeneResults(outputGeneResults, ios::out);
          sgeModel->printResultsForGenes(streamGeneResults);
          streamGeneResults.close();
        }

        // Compute Confidence intervals:
        string outputGeneCIResults = ApplicationTools::getAFilePath("output.results.genes_hessian.file", nisce.getParams(), false, false);
        if (outputGeneCIResults != "none")
        {
          ofstream streamGeneCIResults(outputGeneCIResults, ios::out);
          d1_array_str params;
          params.push_back("r");
          for (auto p : sgeModel->distParamNames_)
          {
            params.push_back(p);
          }
          streamGeneCIResults << "Gene";
          for (auto p : params)
          {
            streamGeneCIResults << "," << p;
          }
          streamGeneCIResults << endl;

          ApplicationTools::displayTask("Computing CI of gene distribution parameters");
          cout << endl;
          progress = 0;
          d2_array_dbl hessianValues(sgeModel->nbGenes_);
          for_each (std::execution::par,
                    sgeModel->geneIndex_.begin(),
                    sgeModel->geneIndex_.end(),
                    [&](size_t i)
          {
            auto f = make_shared<FreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
            auto numDer = make_shared<ThreePointsNumericalDerivative>(f);
            numDer->enableSecondOrderDerivatives(true);
            numDer->enableSecondOrderCrossDerivatives(false);
            numDer->setParametersToDerivate(params);
            numDer->setParameters(numDer->getParameters()); //Induce calculation of derivatives.
            hessianValues[i].resize(params.size());
            for (size_t k = 0; k < params.size(); ++k)
            {
              hessianValues[i][k] = - numDer->getSecondOrderDerivative(params[k]); //We use - here as the function is -logLik
            }
            progress++;
            ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
          });
          ApplicationTools::displayTaskDone();
          for (size_t i = 0; i < sgeModel->nbGenes_; ++i)
          {
            streamGeneCIResults << sgeModel->geneIds_[i];
            for (auto v : hessianValues[i])
            {
              streamGeneCIResults << "," << v;
            }
            streamGeneCIResults << endl; 
          }
          streamGeneCIResults.close();
        }

      }
      else if (method == "LFP" || method == "LessFreeParameters")
      {
        /***********************************************************
        *                                                         *
        * Simple model with free parameters. Gamma or ScaledBeta. *
        * Drop-out parameter has a prior distribution.            *
        *                                                         *
        ***********************************************************/

        // Initialize parameters:
        string model = ApplicationTools::getStringParameter("gene_distribution", nisce.getParams(), "ScaledBeta");
        size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 10);
        size_t nbRCategories = ApplicationTools::getParameter<size_t>("nb_dropout_categories", nisce.getParams(), 10);

        // Create the data structure:
        auto sgeModel = make_shared<LessFreeParameterSgeModel>(counts, geneIds, cellIds, model, nbCategories, nbRCategories);
        sgeModel->print(cout);
        sgeModel->initParameters();

        // Now optimize parameters.
        double tol1 = ApplicationTools::getDoubleParameter("tolerance1", nisce.getParams(), 0.01);

        ApplicationTools::displayTask("Coarse fitting of gene distributions");
        cout << endl;
        unsigned int progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<LessFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol1);
          optimizer->optimize();
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes1 = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
        if (outputGenes1 != "none")
        {
          ofstream streamGenes1(outputGenes1, ios::out);
          sgeModel->printGeneParameters(streamGenes1);
          streamGenes1.close();
        }

        // Hyper parameters: === Note: Beta prior seems difficult to fit, we use a uniform one ===
        //
        //ApplicationTools::displayTask("Coarse fitting of hyper-parameters");
        //cout << endl;
        //auto feb = make_shared<LessFreeParameterSgeModel::FullLikelihood>(sgeModel);
        //auto hyperOptimizer = make_unique<PowellMultiDimensions>(feb);
        //hyperOptimizer->setStopCondition(make_shared<FunctionStopCondition>(hyperOptimizer.get()));
        //hyperOptimizer->setMessageHandler(messenger);
        //hyperOptimizer->setProfiler(profiler);
        //hyperOptimizer->setVerbose(0);
        //hyperOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        //hyperOptimizer->init(feb->getParameters());
        //hyperOptimizer->getStopCondition()->setTolerance(tol1);
        //hyperOptimizer->optimize();
        //ApplicationTools::displayTaskDone();
        //ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());
        //sgeModel->printHyperParameters(cout);  

        // Then cell specific parameters:

        double tol2 = ApplicationTools::getDoubleParameter("tolerance2", nisce.getParams(), 0.0001);

        ApplicationTools::displayTask("Fitting of cell-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->cellIndex_.begin(),
                  sgeModel->cellIndex_.end(),
                  [&](size_t j)
        {
          auto f = make_shared<LessFreeParameterSgeModel::CellLikelihood>(sgeModel, j);
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
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputCells2 = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
        if (outputCells2 != "none")
        {
          ofstream streamCells2(outputCells2, ios::out);
          sgeModel->printCellParameters(streamCells2);
          streamCells2.close();
        }

        //ApplicationTools::displayTask("Re-fitting of hyper-parameters");
        //cout << endl;
        //hyperOptimizer->init(feb->getParameters());
        //hyperOptimizer->getStopCondition()->setTolerance(tol2);
        //hyperOptimizer->optimize();
        //ApplicationTools::displayTaskDone();
        //ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());
        //sgeModel->printHyperParameters(cout);  

        ApplicationTools::displayTask("Re-fitting of gene distributions");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<LessFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes2 = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
        if (outputGenes2 != "none")
        {
          ofstream streamGenes2(outputGenes2, ios::out);
          sgeModel->printGeneParameters(streamGenes2);
          streamGenes2.close();
        }

        // Fine fitting
        
        double tol3 = ApplicationTools::getDoubleParameter("tolerance3", nisce.getParams(), 0.000001);

        ApplicationTools::displayTask("Fine-fitting of cell-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->cellIndex_.begin(),
                  sgeModel->cellIndex_.end(),
                  [&](size_t j) 
        {
          auto f = make_shared<LessFreeParameterSgeModel::CellLikelihood>(sgeModel, j);
          auto optimizer = make_unique<BrentOneDimension>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->setInitialInterval(0.001, 100);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol3);
          optimizer->optimize();
          sgeModel->ls_[j] = optimizer->getParameterValue("l");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputCells3 = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
        if (outputCells3 != "none")
        {
          ofstream streamCells3(outputCells3, ios::out);
          sgeModel->printCellParameters(streamCells3);
          streamCells3.close();
        }

        //ApplicationTools::displayTask("Fine-fitting of hyper-parameters");
        //cout << endl;
        //hyperOptimizer->init(feb->getParameters());
        //hyperOptimizer->getStopCondition()->setTolerance(tol3);
        //hyperOptimizer->optimize();
        //ApplicationTools::displayTaskDone();
        //ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());
        //sgeModel->printHyperParameters(cout);  

        ApplicationTools::displayTask("Fine-fitting of gene distributions");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<LessFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = optimizer->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes3 = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
        if (outputGenes3 != "none")
        {
          ofstream streamGenes3(outputGenes3, ios::out);
          sgeModel->printGeneParameters(streamGenes3);
          streamGenes3.close();
        }

        // Final results:
        string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
        if (outputGeneResults != "none")
        {
          ofstream streamGeneResults(outputGeneResults, ios::out);
          sgeModel->printResultsForGenes(streamGeneResults);
          streamGeneResults.close();
        }

        // Compute Confidence intervals:
        string outputGeneCIResults = ApplicationTools::getAFilePath("output.results.genes_hessian.file", nisce.getParams(), false, false);
        if (outputGeneCIResults != "none")
        {
          ofstream streamGeneCIResults(outputGeneCIResults, ios::out);
          d1_array_str params = sgeModel->distParamNames_;
          streamGeneCIResults << "Gene";
          for (auto p : params)
          {
            streamGeneCIResults << "," << p;
          }
          streamGeneCIResults << endl;

          ApplicationTools::displayTask("Computing CI of gene distribution parameters");
          cout << endl;
          progress = 0;
          d2_array_dbl hessianValues(sgeModel->nbGenes_);
          for_each (std::execution::par,
                    sgeModel->geneIndex_.begin(),
                    sgeModel->geneIndex_.end(),
                    [&](size_t i)
          {
            auto f = make_shared<LessFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
            auto numDer = make_shared<ThreePointsNumericalDerivative>(f);
            numDer->enableSecondOrderDerivatives(true);
            numDer->enableSecondOrderCrossDerivatives(false);
            numDer->setParametersToDerivate(params);
            numDer->setParameters(numDer->getParameters()); //Induce calculation of derivatives.
            hessianValues[i].resize(params.size());
            for (size_t k = 0; k < params.size(); ++k)
            {
              hessianValues[i][k] = - numDer->getSecondOrderDerivative(params[k]); //We use - here as the function is -logLik
            }
            progress++;
            ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
          });
          ApplicationTools::displayTaskDone();
          for (size_t i = 0; i < sgeModel->nbGenes_; ++i)
          {
            streamGeneCIResults << sgeModel->geneIds_[i];
            for (auto v : hessianValues[i])
            {
              streamGeneCIResults << "," << v;
            }
            streamGeneCIResults << endl; 
          }
          streamGeneCIResults.close();
        }

      }

      else if (method == "FEB" || method == "FullEmpiricalBayesian")
      {
        //TODO redo optimization as for the FP model (1,2,3)
        /***********************************************
         *                                             *
         * Full empirical Bayesian:                    *
         * Mean and noise have a-priori distributions. *
         *                                             *
         ***********************************************/

        size_t nbMeanCategories = ApplicationTools::getParameter<size_t>("nb_mean_categories", nisce.getParams(), 5);
        size_t nbShapeCategories = ApplicationTools::getParameter<size_t>("nb_shape_categories", nisce.getParams(), 5);
        size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 5);

        // Create the data structure:
        auto sgeModel = make_shared<FullEmpiricalBayesianSgeModel>(counts, geneIds, cellIds, nbMeanCategories, nbShapeCategories, nbCategories);
        sgeModel->print(cout);
        sgeModel->initParameters();

        // Now optimize parameters.
        double tol1 = ApplicationTools::getDoubleParameter("tolerance1", nisce.getParams(), 0.01);

        // First hyper-parameters:
        ApplicationTools::displayTask("Coarse fitting of hyper-parameters");
        cout << endl;
        auto feb = make_shared<FullEmpiricalBayesianSgeModel::EBLikelihood>(sgeModel);
        auto hyperOptimizer = make_unique<PowellMultiDimensions>(feb);
        hyperOptimizer->setMessageHandler(messenger);
        hyperOptimizer->setProfiler(profiler);
        hyperOptimizer->setVerbose(1);
        hyperOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        hyperOptimizer->init(feb->getParameters());
        hyperOptimizer->getStopCondition()->setTolerance(tol1);
        hyperOptimizer->optimize();
        sgeModel->meanMean_ = hyperOptimizer->getParameterValue("meanMean");
        sgeModel->shapeMean_ = hyperOptimizer->getParameterValue("shapeMean");
        sgeModel->meanShape_ = hyperOptimizer->getParameterValue("meanShape");
        sgeModel->shapeShape_ = hyperOptimizer->getParameterValue("shapeShape");
        ApplicationTools::displayTaskDone();

        sgeModel->printHyperParameters(cout);

        // Then cell specific parameters:
 
        double tol2 = ApplicationTools::getDoubleParameter("tolerance2", nisce.getParams(), 0.0001);

        ApplicationTools::displayTask("Fitting of cell-specific parameters");
        cout << endl;
        unsigned int progress = 0;
        for_each (std::execution::par,
            sgeModel->cellIndex_.begin(),
            sgeModel->cellIndex_.end(),
            [&](size_t j)
       {
          auto f = make_shared<FullEmpiricalBayesianSgeModel::CellLikelihood>(sgeModel, j);
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
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();

        string outputCells1 = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
        if (outputCells1 != "none")
        {
          ofstream streamCells1(outputCells1, ios::out);
          sgeModel->printCellParameters(streamCells1);
          streamCells1.close();
        }

        // Then gene specific parameters:
        ApplicationTools::displayTask("Fitting of gene-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<FullEmpiricalBayesianSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<BrentOneDimension>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->setInitialInterval(0.001, 100);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->rs_[i] = optimizer->getParameterValue("r");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();

        string outputGenes1 = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
        if (outputGenes1 != "none")
        {
          ofstream streamGenes1(outputGenes1, ios::out);
          sgeModel->printGeneParameters(streamGenes1);
          streamGenes1.close();
        }

        // Final hyper-parameters:
        double tol3 = ApplicationTools::getDoubleParameter("tolerance3", nisce.getParams(), 0.000001);
        ApplicationTools::displayTask("Fine fitting of hyper-parameters");
        cout << endl;
        hyperOptimizer->init(feb->getParameters());
        hyperOptimizer->getStopCondition()->setTolerance(tol3);
        hyperOptimizer->optimize();
        sgeModel->meanMean_ = hyperOptimizer->getParameterValue("meanMean");
        sgeModel->shapeMean_ = hyperOptimizer->getParameterValue("shapeMean");
        sgeModel->meanShape_ = hyperOptimizer->getParameterValue("meanShape");
        sgeModel->shapeShape_ = hyperOptimizer->getParameterValue("shapeShape");
        ApplicationTools::displayTaskDone();

        sgeModel->printHyperParameters(cout);

        // Print estimates of mean and noise:
        string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
        if (outputGeneResults != "none")
        {
          ofstream streamGeneResults(outputGeneResults, ios::out);
          sgeModel->printResultsForGenes(streamGeneResults);
          streamGeneResults.close();
        }

        string outputResultMean = ApplicationTools::getAFilePath("output.results.mean.file", nisce.getParams(), false, false);
        ApplicationTools::displayResult("Results for expression mean in", outputResultMean);
        if (outputResultMean != "none")
        {
          ofstream streamResultMean(outputResultMean, ios::out);
          sgeModel->printMeanResults(streamResultMean);
          streamResultMean.close();
        }
        string outputResultNoise = ApplicationTools::getAFilePath("output.results.noise.file", nisce.getParams(), false, false);
        ApplicationTools::displayResult("Results for expression noise in", outputResultNoise);
        if (outputResultNoise != "none")
        {
          ofstream streamResultNoise(outputResultNoise, ios::out);
          sgeModel->printNoiseResults(streamResultNoise);
          streamResultNoise.close();
        }
        // Print entropy of mean and noise:
        string outputResultMeanEntropy = ApplicationTools::getAFilePath("output.results.mean.entropy.file", nisce.getParams(), false, false);
        ApplicationTools::displayResult("Entropy for expression mean in", outputResultMeanEntropy);
        if (outputResultMeanEntropy != "none")
        {
          ofstream streamResultMeanEntropy(outputResultMeanEntropy, ios::out);
          sgeModel->printMeanEntropyResults(streamResultMeanEntropy);
          streamResultMeanEntropy.close();
        }
        string outputResultNoiseEntropy = ApplicationTools::getAFilePath("output.results.noise.entropy.file", nisce.getParams(), false, false);
        ApplicationTools::displayResult("Entropy for expression noise in", outputResultNoiseEntropy);
        if (outputResultNoiseEntropy != "none")
        {
          ofstream streamResultNoiseEntropy(outputResultNoiseEntropy, ios::out);
          sgeModel->printNoiseEntropyResults(streamResultNoiseEntropy);
          streamResultNoiseEntropy.close();
        }
      }
      else if (method == "EB" || method == "EmpiricalBayesian")
      {
        /********************************************
         *                                          *
         * Partial empirical Bayesian:              *
         * Only noise has an a-priori distribution. *
         *                                          *
         ********************************************/

        size_t nbShapeCategories = ApplicationTools::getParameter<size_t>("nb_shape_categories", nisce.getParams(), 5);
        size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 5);

        // Create the data structure:
        auto sgeModel = make_shared<EmpiricalBayesianSgeModel>(counts, geneIds, cellIds, nbShapeCategories, nbCategories);
        sgeModel->print(cout);
        sgeModel->initParameters();

        // Now optimize parameters.
        double tol1 = ApplicationTools::getDoubleParameter("tolerance1", nisce.getParams(), 0.01);

        // First hyper-parameters:
        ApplicationTools::displayTask("Coarse fitting of hyper-parameters");
        cout << endl;
        auto feb = make_shared<EmpiricalBayesianSgeModel::EBLikelihood>(sgeModel);
        auto hyperOptimizer = make_unique<PowellMultiDimensions>(feb);
        hyperOptimizer->setMessageHandler(messenger);
        hyperOptimizer->setProfiler(profiler);
        hyperOptimizer->setVerbose(1);
        hyperOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        hyperOptimizer->init(feb->getParameters());
        hyperOptimizer->getStopCondition()->setTolerance(tol1);
        hyperOptimizer->optimize();
        sgeModel->meanShape_ = hyperOptimizer->getParameterValue("meanShape");
        sgeModel->shapeShape_ = hyperOptimizer->getParameterValue("shapeShape");
        ApplicationTools::displayTaskDone();

        sgeModel->printHyperParameters(cout);

        // Then cell specific parameters:
        vector<size_t> cellIndex;
        for (size_t j = 0; j < sgeModel->nbCells_; ++j)
        {
          cellIndex.push_back(j);
        }

        double tol2 = ApplicationTools::getDoubleParameter("tolerance2", nisce.getParams(), 0.0001);

        ApplicationTools::displayTask("Fitting of cell-specific parameters");
        cout << endl;
        unsigned int progress = 0;
        for_each (std::execution::par, cellIndex.begin(), cellIndex.end(), [&](size_t j) {
          auto f = make_shared<EmpiricalBayesianSgeModel::CellLikelihood>(sgeModel, j);
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
          ApplicationTools::displayGauge(progress, cellIndex.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();

        string outputCells1 = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
        if (outputCells1 != "none")
        {
          ofstream streamCells1(outputCells1, ios::out);
          sgeModel->printCellParameters(streamCells1);
          streamCells1.close();
        }

        // Then gene specific parameters:
        ApplicationTools::displayTask("Fitting of gene-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par, sgeModel->geneIndex_.begin(), sgeModel->geneIndex_.end(), [&](size_t i) {
          auto f = make_shared<EmpiricalBayesianSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->rs_[i] = optimizer->getParameterValue("r");
          sgeModel->means_[i] = optimizer->getParameterValue("mean");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();

        string outputGenes1 = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
        if (outputGenes1 != "none")
        {
          ofstream streamGenes1(outputGenes1, ios::out);
          sgeModel->printGeneParameters(streamGenes1);
          streamGenes1.close();
        }

        // Final hyper-parameters:
        double tol3 = ApplicationTools::getDoubleParameter("tolerance3", nisce.getParams(), 0.000001);
        ApplicationTools::displayTask("Fine fitting of hyper-parameters");
        cout << endl;
        hyperOptimizer->init(feb->getParameters());
        hyperOptimizer->getStopCondition()->setTolerance(tol3);
        hyperOptimizer->optimize();
        sgeModel->meanShape_ = hyperOptimizer->getParameterValue("meanShape");
        sgeModel->shapeShape_ = hyperOptimizer->getParameterValue("shapeShape");
        ApplicationTools::displayTaskDone();

        sgeModel->printHyperParameters(cout);

        // Final gene specific parameters:
        ApplicationTools::displayTask("Fine fitting of gene-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par, sgeModel->geneIndex_.begin(), sgeModel->geneIndex_.end(), [&](size_t i) {
          auto f = make_shared<EmpiricalBayesianSgeModel::GeneLikelihood>(sgeModel, i);
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol3);
          optimizer->optimize();
          sgeModel->rs_[i] = optimizer->getParameterValue("r");
          sgeModel->means_[i] = optimizer->getParameterValue("mean");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();

        string outputGenes2 = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
        if (outputGenes2 != "none")
        {
          ofstream streamGenes2(outputGenes2, ios::out);
          sgeModel->printGeneParameters(streamGenes2);
          streamGenes2.close();
        }

        // Print estimates of mean and noise:
        string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
        if (outputGeneResults != "none")
        {
          ofstream streamGeneResults(outputGeneResults, ios::out);
          sgeModel->printResultsForGenes(streamGeneResults);
          streamGeneResults.close();
        }
      }
    }
    else
    {
      /*********************************************
      *                                            *
      * Read partitions from a file and test them. *******
      * Only the free parameter with Gamma is supported. *
      *                                                  *
      ****************************************************/

      ifstream partitionStream(partitionFile, ios::in);
      ApplicationTools::displayResult("Reading cell partitions from file", partitionFile);
      auto partitions = DataTable::read(partitionStream, TextTools::toString(sep), true, 0);
      size_t nbPartitions = partitions->getNumberOfColumns();
      ApplicationTools::displayResult("Number of partitions found", nbPartitions);

      //Check cells:
      auto cellsInPartitions = partitions->getRowNames();
      for (size_t i = partitions->getNumberOfRows(); i > 0; --i)
      {
        if (!VectorTools::contains(cellIds, cellsInPartitions[i - 1]))
        {
          partitions->deleteRow(i - 1);
          cout << "Cell " << cellsInPartitions[i - 1] << " was ignored." << endl;
        }
      }
      if (partitions->getNumberOfRows() == 0)
      {
        throw Exception("Error!!! The cell names in the partition file do not match those in the count file.");
      }

      // Initialize parameters:
      size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 5);

      double tol1 = ApplicationTools::getDoubleParameter("tolerance1", nisce.getParams(), 0.01);
      double tol2 = ApplicationTools::getDoubleParameter("tolerance2", nisce.getParams(), 0.0001);
      double tol3 = ApplicationTools::getDoubleParameter("tolerance3", nisce.getParams(), 0.000001);
      bool aliasShape = ApplicationTools::getBooleanParameter("alias_shape", nisce.getParams(), true);
      ApplicationTools::displayBooleanResult("Single shape parameter", aliasShape);
      
      string partitionResultFile = ApplicationTools::getStringParameter("cell_partition.result.file", nisce.getParams(), "None");

      if (partitionResultFile == "None" || partitionResultFile == "none")
      {
        ApplicationTools::displayWarning("!!! No partition result will be saved, only estimates. Model comparison will not be possible.");
      }

      // Now evaluate each partition model:
      d1_array_dbl partitionLogLikelihoods(nbPartitions);
      d1_array_dbl partitionNbParameters(nbPartitions);
      for (size_t p = 0; p < nbPartitions; ++p)
      {
        ApplicationTools::displayResult("Fitting partition model", partitions->getColumnName(p));
        // Read partition:
        map<string, std::vector<size_t>> pMap;
        for (size_t j = 0; j < partitions->getNumberOfRows(); ++j) 
        {
          string cell = partitions->getRowName(j);
          pMap[(*partitions)(j, p)].push_back(VectorTools::which(cellIds, cell));
        }
        Partition partition(pMap);

        // Create the data structure:
        auto sgeModel = make_shared<CellPartitionFreeParameterSgeModel>(counts, geneIds, cellIds, partition, nbCategories);
        sgeModel->print(cout);
        sgeModel->initParameters();

        // Now optimize parameters.
        ApplicationTools::displayTask("Coarse fitting of gene distributions");
        cout << endl;
        unsigned int progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                 [&](size_t i)
        {
          auto f = make_shared<CellPartitionFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          if (aliasShape) {
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
          }
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getIndependentParameters());
          optimizer->getStopCondition()->setTolerance(tol1);
          optimizer->optimize();
          sgeModel->rs_[i] = f->getParameterValue("r");
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = f->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes1 = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
        if (outputGenes1 != "none")
        {
          ofstream streamGenes1(outputGenes1 + ".part" + TextTools::toString(p + 1), ios::out);
          sgeModel->printGeneParameters(streamGenes1);
          streamGenes1.close();
        }

        // Then cell specific parameters:

        ApplicationTools::displayTask("Fitting of cell-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->cellIndex_.begin(),
                  sgeModel->cellIndex_.end(),
                  [&](size_t j)
        {
          auto f = make_shared<CellPartitionFreeParameterSgeModel::CellLikelihood>(sgeModel, j);
          auto optimizer = make_unique<BrentOneDimension>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->setInitialInterval(0.001, 100);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->ls_[j] = f->getParameterValue("l");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputCells2 = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
        if (outputCells2 != "none")
        {
          ofstream streamCells2(outputCells2 + ".part" + TextTools::toString(p + 1), ios::out);
          sgeModel->printCellParameters(streamCells2);
          streamCells2.close();
        }

        ApplicationTools::displayTask("Re-fitting of gene distributions");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                 [&](size_t i)
        {
          auto f = make_shared<CellPartitionFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          if (aliasShape) {
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
          }
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getIndependentParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->rs_[i] = f->getParameterValue("r");
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = f->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes2 = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
        if (outputGenes2 != "none")
        {
          ofstream streamGenes2(outputGenes2 + ".part" + TextTools::toString(p + 1), ios::out);
          sgeModel->printGeneParameters(streamGenes2);
          streamGenes2.close();
        }

        // Fine fitting
        
        ApplicationTools::displayTask("Fine-fitting of cell-specific parameters");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->cellIndex_.begin(),
                  sgeModel->cellIndex_.end(),
                  [&](size_t j) 
        {
          auto f = make_shared<CellPartitionFreeParameterSgeModel::CellLikelihood>(sgeModel, j);
          auto optimizer = make_unique<BrentOneDimension>(f);
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->setInitialInterval(0.001, 100);
          optimizer->init(f->getParameters());
          optimizer->getStopCondition()->setTolerance(tol3);
          optimizer->optimize();
          sgeModel->ls_[j] = f->getParameterValue("l");
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->cellIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputCells3 = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
        if (outputCells3 != "none")
        {
          ofstream streamCells3(outputCells3 + ".part" + TextTools::toString(p + 1), ios::out);
          sgeModel->printCellParameters(streamCells3);
          streamCells3.close();
        }

        ApplicationTools::displayTask("Fine-fitting of gene distributions");
        cout << endl;
        progress = 0;
        for_each (std::execution::par,
                  sgeModel->geneIndex_.begin(),
                  sgeModel->geneIndex_.end(),
                  [&](size_t i)
        {
          auto f = make_shared<CellPartitionFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
          if (aliasShape) {
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
          }
          auto optimizer = make_unique<PowellMultiDimensions>(f);
          optimizer->setStopCondition(make_shared<FunctionStopCondition>(optimizer.get()));
          optimizer->setMessageHandler(messenger);
          optimizer->setProfiler(profiler);
          optimizer->setVerbose(0);
          optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
          optimizer->init(f->getIndependentParameters());
          optimizer->getStopCondition()->setTolerance(tol2);
          optimizer->optimize();
          sgeModel->rs_[i] = f->getParameterValue("r");
          for (size_t k = 0; k < sgeModel->nbDistParams_; ++k)
          {
            sgeModel->distParamValues_[i][k] = f->getParameterValue(sgeModel->distParamNames_[k]);
          }
          progress++;
          ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
        });
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayResult("Log likelihood is", sgeModel->getLogLikelihood());

        string outputGenes3 = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
        if (outputGenes3 != "none")
        {
          ofstream streamGenes3(outputGenes3 + ".part" + TextTools::toString(p + 1), ios::out);
          sgeModel->printGeneParameters(streamGenes3);
          streamGenes3.close();
        }

        // Final results:
        string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
        if (outputGeneResults != "none")
        {
          ofstream streamGeneResults(outputGeneResults + ".part" + TextTools::toString(p + 1), ios::out);
          sgeModel->printResultsForGenes(streamGeneResults);
          streamGeneResults.close();
        }

        // Store model results:
        partitionLogLikelihoods[p] = sgeModel->getLogLikelihood();
        partitionNbParameters[p] = sgeModel->getNumberOfParameters();

        // Compute Confidence intervals:
        string outputGeneCIResults = ApplicationTools::getAFilePath("output.results.genes_hessian.file", nisce.getParams(), false, false);
        if (outputGeneCIResults != "none")
        {
          ofstream streamGeneCIResults(outputGeneCIResults + ".part" + TextTools::toString(p + 1), ios::out);
          auto f0 = make_unique<CellPartitionFreeParameterSgeModel::GeneLikelihood>(sgeModel, 0);
          if (aliasShape) {
            auto pnames = f0->getParameters().getParameterNames();
            pnames.erase(remove_if(pnames.begin(),
                                   pnames.end(),
                                   [](const string& name) {
                                     return(!name.ends_with("Shape"));
                                   }),
                         pnames.end());
            for (auto name : pnames) {
              if (name != pnames[0] ) {
                f0->aliasParameters(pnames[0], name);
              }
            }
          }
          d1_array_str params = f0->getIndependentParameters().getParameterNames();
          streamGeneCIResults << "Gene";
          for (auto pm : params)
          {
            streamGeneCIResults << "," << pm;
          }
          streamGeneCIResults << endl;

          ApplicationTools::displayTask("Computing CI of gene distribution parameters");
          cout << endl;
          progress = 0;
          d2_array_dbl hessianValues(sgeModel->nbGenes_);
          for_each (std::execution::par,
                    sgeModel->geneIndex_.begin(),
                    sgeModel->geneIndex_.end(),
                    [&](size_t i)
          {
            auto f = make_shared<CellPartitionFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
            if (aliasShape) {
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
            }
            auto numDer = make_shared<ThreePointsNumericalDerivative>(f);
            numDer->enableSecondOrderDerivatives(true);
            numDer->enableSecondOrderCrossDerivatives(false);
            numDer->setParametersToDerivate(params);
            numDer->setParameters(numDer->getParameters()); //Induce calculation of derivatives.
            hessianValues[i].resize(params.size());
            for (size_t k = 0; k < params.size(); ++k)
            {
              hessianValues[i][k] = - numDer->getSecondOrderDerivative(params[k]); //We use - here as the function is -logLik
            }
            progress++;
            ApplicationTools::displayGauge(progress, sgeModel->geneIndex_.size() - 1, '=');
          });
          ApplicationTools::displayTaskDone();
          for (size_t i = 0; i < sgeModel->nbGenes_; ++i)
          {
            streamGeneCIResults << sgeModel->geneIds_[i];
            for (auto v : hessianValues[i])
            {
              streamGeneCIResults << "," << v;
            }
            streamGeneCIResults << endl; 
          }
          streamGeneCIResults.close();
        }


      }

      // Save model comparison:
      if (partitionResultFile != "None" && partitionResultFile != "none")
      {
        ofstream partitionResultStream(partitionResultFile, ios::out);
        partitionResultStream << "Partition,LogLikelihood,NbParameters,AIC" << endl;
        for (size_t p = 0; p < nbPartitions; ++p)
        {
          partitionResultStream << partitions->getColumnName(p);
          partitionResultStream << "," << partitionLogLikelihoods[p];
          partitionResultStream << "," << partitionNbParameters[p];
          partitionResultStream << "," <<  2 * partitionNbParameters[p] - 2 * partitionLogLikelihoods[p];
          partitionResultStream << endl;
        }
        partitionResultStream.close();
      }
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
