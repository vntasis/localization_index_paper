# Conceptual example analysis

The scripts in this directory were used to generate the results of Figure 1:

- [`simple_example_ratios.awk`](simple_example_ratios.awk): Generates data for
  scenarios like the one in Figure 1A. That is the number of RNA molecules
  coming from the expression of two genes with a specified nucleo-cytosolic
  distribution, and the expected read counts from a theoretic RNAseq
  experiment.
- [`generate_examples.sh`](generate_examples.sh): This script was used to
  produce the data and the barplot of Figure 1C.
