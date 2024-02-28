# ==================================================================
## Convert expression data to json fromat for analysis with cmdstan
# ==================================================================

library(cmdstanr)
load("RData/fpkms_sincseq.RData")

data_dir <- "stan_models/data_sincseq/"

Ngenes <- nrow(fpkms)
data_file <- paste0(data_dir, "sinc_seq.json")

data_list <-
  list(
    N = Ngenes,
    Cyt = fpkms$cytosolic,
    Nuc = fpkms$nuclear,
    Who = fpkms$whole_cell
  )

write_stan_json(data_list, data_file)
