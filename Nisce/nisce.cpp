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
 *        Version:  0.99.0
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
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Exceptions.h>

using namespace std;
using namespace bpp;

#include "ProbFunctions.h"
#include "FreeParameterSgeModel.h"
#include "LessFreeParameterSgeModel.h"
#include "EmpiricalBayesianSgeModel.h"
#include "FullEmpiricalBayesianSgeModel.h"
#include "CellPartitionFreeParameterSgeModel.h"
#include "CellPartitionLessFreeParameterSgeModel.h"


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



void writeCounts(const d2_array_ui& counts, ostream& output,
           const d1_array_str& geneIds, const d1_array_str& cellIds,
           char sep = ',')
{
  output << "Gene";
  for (auto cellId: cellIds) {
    output << sep << cellId;
  }
  output << endl;
  for (size_t i = 0; i < geneIds.size(); ++i)
  {
    output << geneIds[i];    
    for (const auto& it : counts[i])
    {
      output << sep << it;
    }
    output << endl;
  };
}



void help()
{
  cout << "Usage:" << endl;
  cout << "  nisce param=options.bpp" << endl;
}


int main(int argc, char* argv[])
{
  cout << "******************************************************************" << endl;
  cout << "*                  N.I.S.C.E., version 0.99.0                    *" << endl;
  cout << "* Author: J. Dutheil                        Created on  24/04/23 *" << endl;
  cout << "*                                           Last Modif. 30/07/25 *" << endl;
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

    ApplicationTools::displayResult("Reading counts from", inputPath);

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

    string outputSelectedCounts = ApplicationTools::getAFilePath("output.selected_counts.file", nisce.getParams(), false, false);
    if (outputSelectedCounts != "none" && outputSelectedCounts != "None")
    {
      ApplicationTools::displayResult("Output filtered counts to", outputSelectedCounts);
      ofstream outputCounts(outputSelectedCounts, ios::out);
      writeCounts(counts, outputCounts, geneIds, cellIds, sep);
      outputCounts.close();
    }
    
    // Select model:
    string partitionFile = ApplicationTools::getAFilePath("cell_partition.file", nisce.getParams(), false, true);
    string model = ApplicationTools::getStringParameter("model", nisce.getParams(), "FP");

    auto analysis = ApplicationTools::getVectorParameter<string>("analysis", nisce.getParams(), ',', "estimate");
    
    if (partitionFile == "None" || partitionFile == "none")
    {
      // ************************************
      // *                                  *
      // * Fit a single model to all cells. *
      // *                                  *
      // ************************************

      if (model == "FP" || model == "FreeParameters")
      {
        // ***********************************************************
        // *                                                         *
        // * Simple model with free parameters. Gamma or ScaledBeta. *
        // *                                                         *
        // ***********************************************************

        // Initialize parameters:
        string distribution = ApplicationTools::getStringParameter("gene_distribution", nisce.getParams(), "ScaledBeta");
        size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 5);

        // Create the data structure:
        auto sgeModel = make_shared<FreeParameterSgeModel>(counts, geneIds, cellIds, distribution, nbCategories);
        sgeModel->print(cout);
        
        // Init parameter values:
        string initDesc = ApplicationTools::getStringParameter("init_param", nisce.getParams(), "from_counts");
        string initMethod;
        map<string, string> initArgs;
        KeyvalTools::parseProcedure(initDesc, initMethod, initArgs);
        if (initMethod == "from_counts") {
          sgeModel->initParameters();
        } else if (initMethod == "from_files") {
          auto initGenePath = ApplicationTools::getAFilePath("gene_param", initArgs, true, true);
          auto initCellPath = ApplicationTools::getAFilePath("cell_param", initArgs, true, true);
          ApplicationTools::displayResult("Gene init. parameters in", initGenePath);
          ApplicationTools::displayResult("Cell init. parameters in", initCellPath);
          ifstream geneFile(initGenePath, ios::in);
          ifstream cellFile(initCellPath, ios::in);
          auto geneParameters = DataTable::read(geneFile, ",", true, 0);
          auto cellParameters = DataTable::read(cellFile, ",", true, 0);
          sgeModel->initParameters(*geneParameters, *cellParameters);
        } else {
          throw Exception("Unknown init method: " + initMethod);
        }

        if (VectorTools::contains(analysis, "estimate"))
        {
          string prPath = ApplicationTools::getAFilePath("optimization.profiler", nisce.getParams(), false, false);
          string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", nisce.getParams(), false, false);

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

          // Now optimize parameters.

          vector<OptimizationParameters> opt;
          bool optNorm = ApplicationTools::getBooleanParameter("optimization.estimate_norms", nisce.getParams(), true);
          ApplicationTools::displayBooleanResult("Estimate cell norms", optNorm);

          if (optNorm) {
            opt.resize(4);
            opt[0].description = "Initial optimization";
            opt[1].description = "Coarse optimization";
            opt[2].description = "Refining optimization";
            opt[3].description = "Final optimization";
            opt[0].tolerance = 1.0;
            opt[1].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance1", nisce.getParams(), 0.01);
            opt[2].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance2", nisce.getParams(), 0.01);
            opt[3].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance3", nisce.getParams(), 0.01);
            opt[0].estimateNorms = false;
            opt[1].estimateNorms = optNorm;
            opt[2].estimateNorms = optNorm;
            opt[3].estimateNorms = optNorm;
            opt[0].outputCells = "None";
            opt[1].outputCells = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
            opt[2].outputCells = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
            opt[3].outputCells = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
            opt[0].outputGenes = "None";
            opt[1].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
            opt[2].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
            opt[3].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
          } else {
            opt.resize(1);
            opt[0].description = "Optimization";
            opt[0].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", nisce.getParams(), 0.01);
            opt[0].estimateNorms = optNorm;
            opt[0].outputCells = ApplicationTools::getAFilePath("output.parameters.cells.file", nisce.getParams(), false, false);
            opt[0].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes.file", nisce.getParams(), false, false);
          }

          sgeModel -> optimize(profiler, messenger, opt);

          // Final results:
          string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
          if (outputGeneResults != "none" && outputGeneResults != "None")
          {
            ofstream streamGeneResults(outputGeneResults, ios::out);
            sgeModel->printResultsForGenes(streamGeneResults);
            streamGeneResults.close();
          }
        }

        // Compute Confidence intervals:
        if (VectorTools::contains(analysis, "getCI"))
        {
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
            unsigned int progress = 0;
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
        
        //Perform simulations
        if (VectorTools::contains(analysis, "simulate"))
        {
          unsigned int nbSim = ApplicationTools::getParameter<unsigned int>("nb_replicates", nisce.getParams(), 10);
          string outputSimPath = ApplicationTools::getStringParameter("simulation_output_path", nisce.getParams(), "sim%i.csv");
          ApplicationTools::displayResult("Number of replicates", nbSim);
          ApplicationTools::displayTask("Simulating...");
          for (unsigned int k = 0; k < nbSim; ++k) {
            char path[300]; sprintf(path, outputSimPath.c_str(), k);
            ApplicationTools::displayGauge(k, nbSim - 1, '=');
            auto simCounts = sgeModel->simulateCounts();
            //Write to file:
            ofstream outSim(path, ios::out);
            outSim << "Genes";
            for (size_t j = 0; j < cellIds.size(); ++j) {
              outSim << "," << cellIds[j];
            }
            outSim << endl;
            for (size_t i = 0; i < geneIds.size(); ++i) {
              outSim << geneIds[i];
              for (size_t j = 0; j < cellIds.size(); ++j) {
                outSim << "," << (*simCounts)[i][j];
              }
              outSim << endl;
            }
          }
          ApplicationTools::displayTaskDone();
        }

      }
      else if (model == "LFP" || model == "LessFreeParameters")
      {
        // ***********************************************************
        // *                                                         *
        // * Simple model with free parameters. Gamma or ScaledBeta. *
        // * Drop-out parameter has a prior distribution.            *
        // *                                                         *
        // ***********************************************************

        // Initialize parameters:
        string distribution = ApplicationTools::getStringParameter("gene_distribution", nisce.getParams(), "ScaledBeta");
        size_t nbCategories = ApplicationTools::getParameter<size_t>("nb_categories", nisce.getParams(), 10);
        string dropoutPriorDesc = ApplicationTools::getStringParameter("dropout_prior", nisce.getParams(), "Uniform");

        string dropoutPrior;
        map<string, string> dropoutPriorParam;
        KeyvalTools::parseProcedure(dropoutPriorDesc, dropoutPrior, dropoutPriorParam);
        size_t nbRCategories = ApplicationTools::getParameter<size_t>("nb_categories", dropoutPriorParam, 10);

        // Create the data structure:
        auto sgeModel = make_shared<LessFreeParameterSgeModel>(counts, geneIds, cellIds, distribution, nbCategories, nbRCategories, dropoutPrior);
        sgeModel->print(cout);

        // Init parameter values:
        string initDesc = ApplicationTools::getStringParameter("init_param", nisce.getParams(), "from_counts");
        string initMethod;
        map<string, string> initArgs;
        KeyvalTools::parseProcedure(initDesc, initMethod, initArgs);
        if (initMethod == "from_counts") {
          sgeModel->initParameters();
          if (dropoutPrior == "Beta")
          { 
            double a = ApplicationTools::getParameter<double>("alpha", dropoutPriorParam, 0.5);
            double b = ApplicationTools::getParameter<double>("beta", dropoutPriorParam, 1.0);
            sgeModel->setHyperparameters({a, b});
          }
        } else if (initMethod == "from_files") {
          auto initGenePath = ApplicationTools::getAFilePath("gene_param", initArgs, true, true);
          auto initCellPath = ApplicationTools::getAFilePath("cell_param", initArgs, true, true);
          ApplicationTools::displayResult("Gene init. parameters in", initGenePath);
          ApplicationTools::displayResult("Cell init. parameters in", initCellPath);
          ifstream geneFile(initGenePath, ios::in);
          ifstream cellFile(initCellPath, ios::in);
          auto geneParameters = DataTable::read(geneFile, ",", true, 0);
          auto cellParameters = DataTable::read(cellFile, ",", true, 0);
          sgeModel->initParameters(*geneParameters, *cellParameters);
        } else {
          throw Exception("Unknown init method: " + initMethod);
        }
        sgeModel->displayHyperParameters(cout);

        if (VectorTools::contains(analysis, "estimate"))
        { 
          string prPath = ApplicationTools::getAFilePath("optimization.profiler", nisce.getParams(), false, false);
          string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", nisce.getParams(), false, false);

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

          // Now optimize parameters.

          vector<OptimizationParameters> opt;
          bool optNorm = ApplicationTools::getBooleanParameter("optimization.estimate_norms", nisce.getParams(), true);
          bool optHyper = ApplicationTools::getBooleanParameter("optimization.dropout_prior.estimate_hyper_parameters", nisce.getParams(), true);
          ApplicationTools::displayBooleanResult("Estimate cell norms", optNorm);
          ApplicationTools::displayBooleanResult("Estimate hyper-parameters", optHyper);

          if (optNorm || optHyper) {
            opt.resize(3);
            opt[0].description = "Initial optimization";
            opt[1].description = "Refining optimization";
            opt[2].description = "Final optimization";
            opt[0].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance1", nisce.getParams(), 0.01);
            opt[1].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance2", nisce.getParams(), 0.01);
            opt[2].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance3", nisce.getParams(), 0.01);
            opt[0].estimateNorms = optNorm;
            opt[1].estimateNorms = optNorm;
            opt[2].estimateNorms = optNorm;
            opt[0].outputCells = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
            opt[1].outputCells = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
            opt[2].outputCells = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
            opt[0].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
            opt[1].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
            opt[2].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
            opt[0].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes1.file", nisce.getParams(), false, false);
            opt[1].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes2.file", nisce.getParams(), false, false);
            opt[2].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes3.file", nisce.getParams(), false, false);
            opt[0].estimateHyperParameters = optHyper;
            opt[1].estimateHyperParameters = optHyper;
            opt[2].estimateHyperParameters = optHyper;
          } else {
            opt.resize(1);
            opt[0].description = "Optimization";
            opt[0].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", nisce.getParams(), 0.01);
            opt[0].estimateNorms = optNorm;
            opt[0].outputCells = ApplicationTools::getAFilePath("output.parameters.cells.file", nisce.getParams(), false, false);
            opt[0].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes.file", nisce.getParams(), false, false);
            opt[0].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes.file", nisce.getParams(), false, false);
            opt[0].estimateHyperParameters = optHyper;
          }

          sgeModel -> optimize(profiler, messenger, opt);

          // Final results:
          string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
          if (outputGeneResults != "none" && outputGeneResults != "None")
          {
            ofstream streamGeneResults(outputGeneResults, ios::out);
            sgeModel->printResultsForGenes(streamGeneResults);
            streamGeneResults.close();
          }
        }
       
        // Compute Posterior Estimates:
        if (VectorTools::contains(analysis, "getPosteriorEstimates"))
        {
          string outputGenePosteriorResults = ApplicationTools::getAFilePath("output.results.genes_posterior.file", nisce.getParams(), false, false);
          if (outputGenePosteriorResults != "none" && outputGenePosteriorResults != "None")
          {
            ofstream streamGenePosteriorResults(outputGenePosteriorResults, ios::out);
            sgeModel->printGenePosteriorEstimates(streamGenePosteriorResults);
            streamGenePosteriorResults.close();
          }
        }

        // Compute Confidence intervals:
        if (VectorTools::contains(analysis, "getCI"))
        {
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
            unsigned int progress = 0;
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
      } else {
        throw Exception("Unsupported model: " + model);
      }
    }
    else
    {
      // *********************************************
      // *                                            *
      // * Read partitions from a file and test them. *
      // *                                            *
      // **********************************************

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

      bool aliasShape = ApplicationTools::getBooleanParameter("alias_shape", nisce.getParams(), false);
      ApplicationTools::displayBooleanResult("Single shape parameter", aliasShape);
      if (aliasShape) {
        ApplicationTools::displayWarning("Shape parameters are aliased, degrees of freedom and AIC values will be incorrect in output files.");
      }
      
      string partitionResultFile = ApplicationTools::getStringParameter("cell_partition.result.file", nisce.getParams(), "None");

      if (partitionResultFile == "None" || partitionResultFile == "none")
      {
        ApplicationTools::displayWarning("!!! No partition result will be saved, only estimates. Model comparison will not be possible.");
      }

      string prPath = ApplicationTools::getAFilePath("optimization.profiler", nisce.getParams(), false, false);
      string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", nisce.getParams(), false, false);

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

        if (model == "FP" || model == "FreeParameters")
        {
          // ***********************************************************
          // *                                                         *
          // * Simple model with free parameters. Gamma or ScaledBeta. *
          // *                                                         *
          // ***********************************************************
          
          // Create the data structure:
          auto sgeModel = make_shared<CellPartitionFreeParameterSgeModel>(counts, geneIds, cellIds, partition, nbCategories);
          sgeModel->print(cout);

          // Init parameter values:
          string initDesc = ApplicationTools::getStringParameter("init_param", nisce.getParams(), "from_counts");
          string initMethod;
          map<string, string> initArgs;
          KeyvalTools::parseProcedure(initDesc, initMethod, initArgs);
          if (initMethod == "from_counts") {
            sgeModel->initParameters();
          } else if (initMethod == "from_files") {
            throw Exception("Not yet implemented: " + initMethod);
            //auto initGenePath = ApplicationTools::getAFilePath("gene_param", initArgs, true, true);
            //auto initCellPath = ApplicationTools::getAFilePath("cell_param", initArgs, true, true);
            //ApplicationTools::displayResult("Gene init. parameters in", initGenePath);
            //ApplicationTools::displayResult("Cell init. parameters in", initCellPath);
            //ifstream geneFile(initGenePath, ios::in);
            //ifstream cellFile(initCellPath, ios::in);
            //auto geneParameters = DataTable::read(geneFile, ",", true, 0);
            //auto cellParameters = DataTable::read(cellFile, ",", true, 0);
            //sgeModel->initParameters(*geneParameters, *cellParameters);
          } else {
            throw Exception("Unknown init method: " + initMethod);
          }
        
          if (VectorTools::contains(analysis, "estimate"))
          { 
            vector<OptimizationParameters> opt;
            bool optNorm = ApplicationTools::getBooleanParameter("optimization.estimate_norms", nisce.getParams(), true);
            ApplicationTools::displayBooleanResult("Estimate cell norms", optNorm);

            if (optNorm) {
              opt.resize(4);
              opt[0].description = "Initial optimization";
              opt[1].description = "Coarse optimization";
              opt[2].description = "Refining optimization";
              opt[3].description = "Final optimization";
              opt[0].tolerance = 1.0;
              opt[1].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance1", nisce.getParams(), 0.01);
              opt[2].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance2", nisce.getParams(), 0.01);
              opt[3].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance3", nisce.getParams(), 0.01);
              opt[0].estimateNorms = false;
              opt[1].estimateNorms = optNorm;
              opt[2].estimateNorms = optNorm;
              opt[3].estimateNorms = optNorm;
              opt[0].outputCells = "None";
              opt[1].outputCells = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
              if (opt[1].outputCells != "none" && opt[1].outputCells != "None")
                opt[1].outputCells += ".part" + TextTools::toString(p + 1);
              opt[2].outputCells = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
              if (opt[2].outputCells != "none" && opt[2].outputCells != "None")
                opt[2].outputCells += ".part" + TextTools::toString(p + 1);
              opt[3].outputCells = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
              if (opt[3].outputCells != "none" && opt[3].outputCells != "None")
                opt[3].outputCells += ".part" + TextTools::toString(p + 1);
              opt[0].outputGenes = "None";
              opt[1].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
              if (opt[1].outputGenes != "none" && opt[1].outputGenes != "None")
                opt[1].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[2].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
              if (opt[2].outputGenes != "none" && opt[2].outputGenes != "None")
                opt[2].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[3].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
              if (opt[3].outputGenes != "none" && opt[3].outputGenes != "None")
                opt[3].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[0].aliasShape = aliasShape;
              opt[1].aliasShape = aliasShape;
              opt[2].aliasShape = aliasShape;
              opt[3].aliasShape = aliasShape;
            } else {
              opt.resize(1);
              opt[0].description = "Optimization";
              opt[0].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", nisce.getParams(), 0.01);
              opt[0].estimateNorms = optNorm;
              opt[0].outputCells = ApplicationTools::getAFilePath("output.parameters.cells.file", nisce.getParams(), false, false);
              if (opt[0].outputCells != "none" && opt[0].outputCells != "None")
                opt[0].outputCells += ".part" + TextTools::toString(p + 1);
              opt[0].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes.file", nisce.getParams(), false, false);
              if (opt[0].outputGenes != "none" && opt[0].outputGenes != "None")
              opt[0].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[0].aliasShape = aliasShape;
            }

            sgeModel -> optimize(profiler, messenger, opt);

            // Final results:
            string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
            if (outputGeneResults != "none" && outputGeneResults != "None")
            {
              ofstream streamGeneResults(outputGeneResults + ".part" + TextTools::toString(p + 1), ios::out);
              sgeModel->printResultsForGenes(streamGeneResults);
              streamGeneResults.close();
            }

            // Store model results:
            partitionLogLikelihoods[p] = sgeModel->getLogLikelihood();
            partitionNbParameters[p] = sgeModel->getNumberOfParameters();
	  }

          if (VectorTools::contains(analysis, "getCI"))
          { 
            // Compute Confidence intervals:
            string outputGeneCIResults = ApplicationTools::getAFilePath("output.results.genes_hessian.file", nisce.getParams(), false, false);
            if (outputGeneCIResults != "none" && outputGeneCIResults != "None")
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
              unsigned int progress = 0;
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

        } else if (model == "LFP" || model == "LessFreeParameters") {

          // ***********************************************************
          // *                                                         *
          // * Simple model with free parameters. Gamma or ScaledBeta. *
          // * Drop-out parameter has a prior distribution.            *
          // *                                                         *
          // ***********************************************************

          // Initialize parameters:
          string dropoutPriorDesc = ApplicationTools::getStringParameter("dropout_prior", nisce.getParams(), "Uniform");
          string dropoutPrior;
          map<string, string> dropoutPriorParam;
          KeyvalTools::parseProcedure(dropoutPriorDesc, dropoutPrior, dropoutPriorParam);
          size_t nbRCategories = ApplicationTools::getParameter<size_t>("nb_categories", dropoutPriorParam, 10);
         
          // Create the data structure:
          auto sgeModel = make_shared<CellPartitionLessFreeParameterSgeModel>(counts, geneIds, cellIds, partition, nbCategories, nbRCategories, dropoutPrior);
          sgeModel->print(cout);

          // Init parameter values:
          string initDesc = ApplicationTools::getStringParameter("init_param", nisce.getParams(), "from_counts");
          string initMethod;
          map<string, string> initArgs;
          KeyvalTools::parseProcedure(initDesc, initMethod, initArgs);
          if (initMethod == "from_counts") {
            sgeModel->initParameters();
            if (dropoutPrior == "Beta")
            { 
              double a = ApplicationTools::getParameter<double>("alpha", dropoutPriorParam, 0.5);
              double b = ApplicationTools::getParameter<double>("beta", dropoutPriorParam, 1.0);
              sgeModel->setHyperparameters({a, b});
            }
          } else if (initMethod == "from_files") {
            throw Exception("Not yet implemented: " + initMethod);
            //auto initGenePath = ApplicationTools::getAFilePath("gene_param", initArgs, true, true);
            //auto initCellPath = ApplicationTools::getAFilePath("cell_param", initArgs, true, true);
            //ApplicationTools::displayResult("Gene init. parameters in", initGenePath);
            //ApplicationTools::displayResult("Cell init. parameters in", initCellPath);
            //ifstream geneFile(initGenePath, ios::in);
            //ifstream cellFile(initCellPath, ios::in);
            //auto geneParameters = DataTable::read(geneFile, ",", true, 0);
            //auto cellParameters = DataTable::read(cellFile, ",", true, 0);
            //sgeModel->initParameters(*geneParameters, *cellParameters);
          } else {
            throw Exception("Unknown init method: " + initMethod);
          }
          sgeModel->displayHyperParameters(cout);

          if (VectorTools::contains(analysis, "estimate"))
          { 
            vector<OptimizationParameters> opt;
            bool optNorm = ApplicationTools::getBooleanParameter("optimization.estimate_norms", nisce.getParams(), true);
            bool optHyper = ApplicationTools::getBooleanParameter("optimization.dropout_prior.estimate_hyper_parameters", nisce.getParams(), true);
            ApplicationTools::displayBooleanResult("Estimate cell norms", optNorm);
            ApplicationTools::displayBooleanResult("Estimate hyper-parameters", optHyper);
    
            if (optNorm || optHyper) {
              opt.resize(3);
              opt[0].description = "Initial optimization";
              opt[1].description = "Refining optimization";
              opt[2].description = "Final optimization";
              opt[0].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance1", nisce.getParams(), 0.01);
              opt[1].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance2", nisce.getParams(), 0.01);
              opt[2].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance3", nisce.getParams(), 0.01);
              opt[0].estimateNorms = optNorm;
              opt[1].estimateNorms = optNorm;
              opt[2].estimateNorms = optNorm;
              opt[0].outputCells = ApplicationTools::getAFilePath("output.parameters.cells1.file", nisce.getParams(), false, false);
              if (opt[0].outputCells != "none" && opt[0].outputCells != "None")
                opt[0].outputCells += ".part" + TextTools::toString(p + 1);
               opt[1].outputCells = ApplicationTools::getAFilePath("output.parameters.cells2.file", nisce.getParams(), false, false);
              if (opt[1].outputCells != "none" && opt[1].outputCells != "None")
                opt[1].outputCells += ".part" + TextTools::toString(p + 1);
              opt[2].outputCells = ApplicationTools::getAFilePath("output.parameters.cells3.file", nisce.getParams(), false, false);
              if (opt[2].outputCells != "none" && opt[2].outputCells != "None")
                opt[2].outputCells += ".part" + TextTools::toString(p + 1);
              opt[0].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes1.file", nisce.getParams(), false, false);
              if (opt[0].outputGenes != "none" && opt[0].outputGenes != "None")
                opt[0].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[1].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes2.file", nisce.getParams(), false, false);
              if (opt[1].outputGenes != "none" && opt[1].outputGenes != "None")
                opt[1].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[2].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes3.file", nisce.getParams(), false, false);
              if (opt[2].outputGenes != "none" && opt[2].outputGenes != "None")
                opt[2].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[0].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes1.file", nisce.getParams(), false, false);
              if (opt[0].outputGenesHyper != "none" && opt[0].outputGenesHyper != "None")
                opt[0].outputGenesHyper += ".part" + TextTools::toString(p + 1);
              opt[1].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes2.file", nisce.getParams(), false, false);
              if (opt[1].outputGenesHyper != "none" && opt[1].outputGenesHyper != "None")
                opt[1].outputGenesHyper += ".part" + TextTools::toString(p + 1);
              opt[2].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes3.file", nisce.getParams(), false, false);
              if (opt[2].outputGenesHyper != "none" && opt[2].outputGenesHyper != "None")
                opt[2].outputGenesHyper += ".part" + TextTools::toString(p + 1);
              opt[0].estimateHyperParameters = optHyper;
              opt[1].estimateHyperParameters = optHyper;
              opt[2].estimateHyperParameters = optHyper;
              opt[0].aliasShape = aliasShape;
              opt[1].aliasShape = aliasShape;
              opt[2].aliasShape = aliasShape;
            } else {
              opt.resize(1);
              opt[0].description = "Optimization";
              opt[0].tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", nisce.getParams(), 0.01);
              opt[0].estimateNorms = optNorm;
              opt[0].outputCells = ApplicationTools::getAFilePath("output.parameters.cells.file", nisce.getParams(), false, false);
              if (opt[0].outputCells != "none" && opt[0].outputCells != "None")
                opt[0].outputCells += ".part" + TextTools::toString(p + 1);
              opt[0].outputGenes = ApplicationTools::getAFilePath("output.parameters.genes.file", nisce.getParams(), false, false);
              if (opt[0].outputGenes != "none" && opt[0].outputGenes != "None")
                opt[0].outputGenes += ".part" + TextTools::toString(p + 1);
              opt[0].outputGenesHyper = ApplicationTools::getAFilePath("output.hyper-parameters.genes.file", nisce.getParams(), false, false);
              if (opt[0].outputGenesHyper != "none" && opt[0].outputGenesHyper != "None")
                opt[0].outputGenesHyper += ".part" + TextTools::toString(p + 1);
              opt[0].estimateHyperParameters = optHyper;
              opt[0].aliasShape = aliasShape;
            }

            sgeModel -> optimize(profiler, messenger, opt);

            // Final results:
            string outputGeneResults = ApplicationTools::getAFilePath("output.results.genes.file", nisce.getParams(), false, false);
            if (outputGeneResults != "none" && outputGeneResults != "None")
            {
              ofstream streamGeneResults(outputGeneResults + ".part" + TextTools::toString(p + 1), ios::out);
              sgeModel->printResultsForGenes(streamGeneResults);
              streamGeneResults.close();
            }
            
	    // Store model results:
            partitionLogLikelihoods[p] = sgeModel->getLogLikelihood();
            partitionNbParameters[p] = sgeModel->getNumberOfParameters();
	  }

          // Compute Posterior Estimates:
          if (VectorTools::contains(analysis, "getPosteriorEstimates"))
          {
            string outputGenePosteriorResults = ApplicationTools::getAFilePath("output.results.genes_posterior.file", nisce.getParams(), false, false);
            if (outputGenePosteriorResults != "none" && outputGenePosteriorResults != "None")
            {
              ofstream streamGenePosteriorResults(outputGenePosteriorResults, ios::out);
              sgeModel->printGenePosteriorEstimates(streamGenePosteriorResults);
              streamGenePosteriorResults.close();
            }
          }

          if (VectorTools::contains(analysis, "getCI"))
          { 
            // Compute Confidence intervals:
            string outputGeneCIResults = ApplicationTools::getAFilePath("output.results.genes_hessian.file", nisce.getParams(), false, false);
            if (outputGeneCIResults != "none")
            {
              ofstream streamGeneCIResults(outputGeneCIResults + ".part" + TextTools::toString(p + 1), ios::out);
              auto f0 = make_unique<CellPartitionLessFreeParameterSgeModel::GeneLikelihood>(sgeModel, 0);
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
              unsigned int progress = 0;
              d2_array_dbl hessianValues(sgeModel->nbGenes_);
              for_each (std::execution::par,
                        sgeModel->geneIndex_.begin(),
                        sgeModel->geneIndex_.end(),
                        [&](size_t i)
              {
                auto f = make_shared<CellPartitionLessFreeParameterSgeModel::GeneLikelihood>(sgeModel, i);
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

        } else {
          throw Exception("Unsupported model: " + model);
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
