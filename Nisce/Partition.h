// SPDX-FileCopyrightText: 2024 Julien Y. Dutheil
//
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * =====================================================================================
 *
 *       Filename:  Partition.h
 *
 * =====================================================================================
 */

#ifndef _PARTITION_H_
#define _PARTITION_H_

class Partition
{
public:

  std::map<std::string, std::vector<size_t>> partitionMap; //Partition name with corresponding cell indexes.

  Partition(const std::map<std::string, std::vector<size_t>>& pMap):
    partitionMap(pMap)
  {}

  size_t size() const
  {
    return partitionMap.size();
  }

  size_t getPartitionIndex(size_t cell) const
  {
    size_t partitionIndex = 0;
    for (auto& p : partitionMap)
    {
      if (std::find(p.second.begin(), p.second.end(), cell) != p.second.end())
      {
        return partitionIndex;
      }
      partitionIndex++;
    }
    throw Exception("Partition::getPartitionIndex(cell): given cell is not part of the specified partition.");
  }

  const std::vector<size_t>& getPartition(size_t partitionIndex) const
  {
    auto it = partitionMap.begin();
    while (partitionIndex > 0) {
      it++;
      partitionIndex--;
    }
    return it->second;
  }
};

#endif //_PARTITION_H_

