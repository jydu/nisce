<!--
SPDX-FileCopyrightText: 2024 Julien Y. Dutheil

SPDX-License-Identifier: GPL-3.0-or-later
-->


# Example usage

The input data is in file `S_paradoxus_CBS8841_rep1_filtered_nisce.csv.gz`.
The file contains one gene per line and one cell per column. Entries are separated by comas, but other characters can be specified.
Input files can be gzipped. The first column contains gene names, and the first line contains cell names.

## Simple model

```bash
nisce param=nisce.bpp
```

Output files:

- `S_paradoxus_CBS8841_rep1_Gamma_cell_parameters.cs`: estimated cell parameters (currently only normalization constants).
- `S_paradoxus_CBS8841_rep1_Gamma_gene_parameters.csv`: gene parameters, also including profile likelihood and profile AIC. 
- `S_paradoxus_CBS8841_rep1_Gamma_gene_results.csv`: gene results (a more user-friendly version of the gene parameters, notably containing noise estimates).
- `S_paradoxus_CBS8841_rep1_Gamma_gene_hessian.csv`: optional, the second order derivatives for each gene parameter. Can be used to get confidence intervals.
- `S_paradoxus_CBS8841_rep1_Gamma_optimization_messages.log`: optional, optimization debugging infos
- `S_paradoxus_CBS8841_rep1_Gamma_optimization_profile.log`: optional, optimization debugging infos


## Partition models

The partition model requires a second input file, listing one or more partitions (`Clusters_S_paradoxus_CBS8841_rep1_filtered_Seurat.csv`). 
The partition file is a simple CSV file with one line per cell, and one column per partititon scheme.
The first column contains cell names, which should match the input count file.
The first line contains names of each partition scheme. 
Partitions are integers: cells with the same number are in the same cluster.

```bash
nisce param=nisce_partition.bpp
```

The output file are the same as for the simple model, excepted that there are cell and gene parameter files for each partition scheme (with file extensions `.part1`, `.part2`, etc.
There is also a file named `S_paradoxus_CBS8841_rep1_Gamma_partition_results.csv`, which contains summary infos for each partition model (including AIC values).
