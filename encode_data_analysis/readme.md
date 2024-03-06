# ENCODE data analysis

Analysis of the ENCODE whole-cell, nuclear, cytosolic RNAseq data

- [main.R](main.R) is the main script of the analysis. It includes
  code to generate Figures 3A and 3B.
- [apply_npidr.R](apply_npidr.R) applies the `npidr` function on the data.
- [npidr.R](npidr.R) defines a function for the non-parametric Irreproducible
  Discovery Rate.
- [restructure_counts.R](restructure_counts.R) is called by the main
  script and helps reorganize the simulated expression data table into a list.
- [preprocess.R](preprocess.R) is called by the main script and
  helps with the filtering and normalization of the data.
- [convert2json.R](convert2json.R) saves the fpkm expression values
  in a json format. This is the input for the
  [inference of _beta_](../beta_inference/) using Stan-NF.
- [save_estimations.R](save_estimations.R) reads the output of
  Stan-NF and saves the posterior estimates in a `.RData` file.
- [feature_association.R](feature_association.R) calls the scripts related to
  the association of the localization index with different transcript features,
  results regarding with the plots for Figure 4.
