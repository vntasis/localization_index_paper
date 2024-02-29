#*****************************
## Filtering and Normalization
#*****************************

#-------------------------------------------------------
## Filter for npIDR < 0.1
## in all out of
## nuc, cyt, who
## Filter low cpm transcripts
## Requirement: The whole cell fraction
## and at least one of the other two (nuclear, cytosolic)
## to have cpm greater than 1
#--------------------------------------------------------


for (sample in names(counts_list)) {
  CPM <-
    counts_list[[sample]][, -1] %>%
    cpm()

  filter_1 <-
    CPM %>%
    `>=`(1) %>%
    rowSums() %>%
    `>`(1)

  filter_2 <-
    CPM[, "whole_cell"] %>%
    `>=`(1)

  if (sample != "H1_2") {
    filter <- filter_1 & filter_2 & np_idr_score_filter[[sample]]
  } else {
    filter <- filter_1 & filter_2
  }


  counts_list[[sample]] %<>% slice(which(filter)) %>%
    slice(grep("Spike", .$transcript_id, invert = TRUE))
  lengths_list[[sample]] %<>% slice(which(filter)) %>%
    slice(grep("Spike", .$transcript_id, invert = TRUE))
}



#--------------------------------------
## Filter mitochondrial transcripts out
#--------------------------------------

for (sample in names(counts_list)) {
  idxs <-
    counts_list[[sample]] %$%
    transcript_id %>%
    match(annotation$transcript_id)

  counts_list[[sample]] %<>%
    add_column(
      .before = 1,
      chr = annotation$chr[idxs]
    ) %>%
    filter(!(chr %in% c("chrM"))) %>%
    select(-chr)

  lengths_list[[sample]] %<>%
    add_column(
      .before = 1,
      chr = annotation$chr[idxs]
    ) %>%
    filter(!(chr %in% c("chrM"))) %>%
    select(-chr)
}



#------------------
## Calculate FPKMs
#------------------

fpkms <- list()
for (sample in names(counts_list)) {

  lib_size <-
    counts_list[[sample]][, -1] %>%
    colSums()

  fpkms[[sample]] <-
    `/`(counts_list[[sample]][, -1], lengths_list[[sample]][, -1]) %>%
    as.matrix() %>%
    `*`(1e+9) %>%
    `%*%`(diag(1 / lib_size)) %>%
    as_tibble() %>%
    add_column(
      .before = 1,
      transcript_id = counts_list[[sample]]$transcript_id
    )

  colnames(fpkms[[sample]]) <- colnames(counts_list[[sample]])

  # discard NA
  idxs <-
    fpkms[[sample]] %>%
    is.na() %>%
    apply(2, which) %>%
    unlist(use.names = FALSE) %>%
    unique()

  if (length(idxs) > 0) fpkms[[sample]] %<>% slice(-idxs)
}


#----------
## Clean up
#----------

rm(list = c("filter_1", "filter_2", "filter",
            "idxs", "sample", "lib_size", "CPM"))
