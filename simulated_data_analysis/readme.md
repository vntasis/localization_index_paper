# Simulated Data Analysis

Scripts used in the generation and analysis of the simulated data:

- [`simulate_whole_cell.sh`](simulation/scripts/simulate_whole_cell.sh) and
  [`simulate_nuclear_cytosolic.sh`](simulation/scripts/simulate_nuclear_cytosolic.sh)
  scripts were used to run the
  [Flux Simulator](https://confluence.sammeth.net/display/SIM/Home).
- The [flux_parameters](simulation/flux_parameters/) directory contains the
  parameter files used with Flux for the generation of whole-cell, nuclear, and
  cytosolic simulated RNAseq data.
- [split2sumcomp.R](simulation/scripts/split2sumcomp.R) was used to split the
  number of molecules per transcript of whole-cell data into nuclear and
  cytosolic.
- [main.R](processing/main.R) is the main script of the analysis. It includes
  code to generate the plots of Figure 2.
- [restructure_counts.R](processing/restructure_counts.R) is called by the main
  script and helps reorganize the simulated expression data table into a list.
- [preprocess.R](processing/preprocess.R) is called by the main script and
  helps with the filtering and normalization of the simulated expression data.
- [convert2json.R](processing/convert2json.R) saves the fpkm expression values
  in a json format. This is the input for the
  [inference of _beta_](../beta_inference/) using Stan-NF.
- [save_estimations.R](processing/save_estimations.R) reads the output of
  Stan-NF and saves the posterior estimates in a `.RData` file.
