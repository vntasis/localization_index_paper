#-------------------------
## Load required libraries
#-------------------------

library(rstan)
library(bayestestR)

load("RData/data_sincseq.RData")

#------------------------------------------------------------
## Fetch posterior estimations of the cytoplasmic proportions
#------------------------------------------------------------

chains <- 4

paths <- paste0(
  "stan_models/results_sincseq/sinc_seq/samples/sinc_seq",
  as.character(1:chains), ".csv"
)

fit <- read_stan_csv(paths)
fit_array <- as.array(fit)
betas <- fit_array[, , "beta"]

posterior <- as.numeric(map_estimate(betas))



#----------------------------------------------
## Save posterior estimations in an .RData file
#----------------------------------------------
save.image(file = "RData/posterior_estimations_sincseq_data.RData")
