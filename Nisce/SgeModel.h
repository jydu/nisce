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
    ls_(nbCells_)
  {}

};

#endif //_SGE_MODEL_H_

