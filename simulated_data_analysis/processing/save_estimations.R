#-------------------------
## Load required libraries
#-------------------------

library(rstan)
library(bayestestR)

load("RData/simulated_data_final.RData")

#------------------------------------------------------------
## Fetch posterior estimations of the cytoplasmic proportions
#------------------------------------------------------------

samples <- names(fpkms)

chains <- 4

posterior <- list()
for (sample in samples) {
  paths <- paste0(
    "stan_models/results_simulation_data/", sample, "/samples/", sample,
    as.character(1:chains), ".csv"
  )

  fit <- read_stan_csv(paths)
  fit_array <- as.array(fit)
  betas <- fit_array[, , "beta"]

  posterior[[sample]][["beta"]] <- as.numeric(map_estimate(betas))
}



#----------------------------------------------
## Save posterior estimations in an .RData file
#----------------------------------------------
save(list = c("posterior"),
     file = "RData/posterior_estimations_simulation_data.RData")
