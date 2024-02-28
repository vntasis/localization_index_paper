# ==================================================================
## Convert expression data to json fromat for analysis with cmdstan
# ==================================================================

library(cmdstanr)
load("RData/simulated_data_final.RData")

data_dir <- "stan_models/data_simulated/"

for (sample in names(fpkms)) {
  data_table <- fpkms[[sample]]
  Ngenes <- nrow(data_table)
  data_file <- paste0(data_dir, sample, ".json")

  data_list <-
    list(
      N = Ngenes,
      Cyt = data_table$cytosolic,
      Nuc = data_table$nuclear,
      Who = data_table$whole_cell
    )

  write_stan_json(data_list, data_file)
}
