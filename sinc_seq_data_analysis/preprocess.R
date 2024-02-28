#*****************************
## Filtering and Normalization
#*****************************

#-------------------------------------------------------
## Filter low cpm genes
## Requirement: The whole cell fraction
## and at least one of the other two (nuclear, cytosolic)
## to have cpm greater than 1
#--------------------------------------------------------

CPM <-
  counts[, -1] %>%
  cpm()


filter_1 <-
  CPM %>%
  `>=`(1) %>%
  rowSums() %>%
  `>`(1)

filter_2 <-
  CPM[, "whole_cell"] %>%
  `>=`(1)

filter <- filter_1 & filter_2

counts %<>% slice(which(filter))
lengths %<>% slice(which(filter))


#--------------------------------------
## Filter mitochondrial genes out
#--------------------------------------

idxs <-
  counts %$%
  transcript_id %>%
  match(annotation$transcript_id)

counts %<>%
  add_column(
    .before = 1,
    chr = annotation$chr[idxs]
  ) %>%
  filter(!(chr %in% c("chrM"))) %>%
  select(-chr)

lengths %<>%
  add_column(
    .before = 1,
    chr = annotation$chr[idxs]
  ) %>%
  filter(!(chr %in% c("chrM"))) %>%
  select(-chr)



#------------------
## Calculate FPKMs
#------------------

lib_size <-
  counts[, -1] %>%
  colSums()

fpkms <-
  `/`(counts[, -1], lengths[, -1]) %>%
  as.matrix() %>%
  `*`(1e+9) %>%
  `%*%`(diag(1 / lib_size)) %>%
  as_tibble() %>%
  add_column(
    .before = 1,
    transcript_id = counts$transcript_id
  )

colnames(fpkms) <- colnames(counts)

# discard NA
idxs <-
  fpkms %>%
  is.na() %>%
  apply(2, which) %>%
  unlist(use.names = FALSE) %>%
  unique()


if (length(idxs) > 0) fpkms %<>% slice(-idxs)


#----------
## Clean up
#----------

rm(list = c("filter_1", "filter_2", "filter", "idxs", "lib_size", "CPM"))
