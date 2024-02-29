#***************************************************
# Apply npIDR to discard irreproducible transcripts*
#***************************************************

source("scripts/npidr.R")


cell_types <- meta %>%
  pull("biosample") %>%
  unique()
fractions <- meta %>%
  pull("fraction") %>%
  unique()

np_idr_score <- list()

for (cell_type in cell_types) {
  for (frac in fractions) {
    samples <-
      meta %>%
      filter(biosample == cell_type & fraction == frac) %>%
      pull(file)

    if (length(samples) > 2) {
      samples <-
        meta %>%
        filter(biosample == cell_type & fraction == frac) %>%
        filter(replicate != 5) %>%
        pull(file)
    }
    if (length(samples) < 2) next

    idr <-
      counts %>%
      select(all_of(samples)) %>%
      npidr(10)

    for (s in samples) np_idr_score[[s]] <- idr
  }
}

# Clean
rm(list = c("idr", "s", "frac", "cell_type",
            "samples", "cell_types", "fractions"))
