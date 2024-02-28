#*****************************
## Filtering and Normalization
#*****************************

fpkms <- list()

for (split in splits) {
  ## Filter low cpm transcripts
  counts_split <- counts_list[[split]]
  lengths_split <- lengths_list[[split]]

  filter_1 <-
    counts_split %>%
    column_to_rownames("transcript_id") %>%
    cpm() %>%
    `>=`(1) %>%
    rowSums() %>%
    `>`(1)

  filter_2 <-
    counts_split %>%
    column_to_rownames("transcript_id") %>%
    cpm() %>%
    `[`(, "whole_cell") %>%
    `>=`(1)

  filter <- filter_1 & filter_2

  counts_split %<>% slice(which(filter))
  lengths_split %<>%
    slice(match(counts_split$transcript_id, lengths_split$transcript_id))



  ## Calculate FPKMs
  lib_size <-
    counts_split[, -1] %>%
    colSums()

  fpkms[[split]] <-
    `/`(counts_split[, -1], lengths_split[, -1]) %>%
    as.matrix() %>%
    `*`(1e+9) %>%
    `%*%`(diag(1 / lib_size)) %>%
    as_tibble() %>%
    add_column(
      .before = 1,
      transcript_id = counts_split$transcript_id
    )

  colnames(fpkms[[split]]) <- colnames(counts_split)

  # discard NA
  idxs <-
    fpkms[[split]] %>%
    is.na() %>%
    apply(2, which) %>%
    unlist(use.names = FALSE) %>%
    unique()

  if (length(idxs) > 0) fpkms[[split]] %<>% slice(-idxs)
}

# Clean up
rm(list = c(
  "filter_1", "filter_2", "filter",
  "counts_split", "lengths_split",
  "split", "idxs", "lib_size"
))
