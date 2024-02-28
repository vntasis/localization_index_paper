#*****************************************************************************
## Turn the count, lenths, and num.of.molecules tibbles into a list of tibble.
## Each element of the list will correspond to
## whole cell, cytosolic, nuclear counts of a sample.
#*****************************************************************************

counts_list <- list()
lengths_list <- list()
n_molecules_list <- list()
for (split in splits) {
  counts_list[[split]] <-
    counts %>%
    select(
      transcript_id, whole_cell,
      paste0(split, "_nuclear"), paste0(split, "_cytosolic")
    ) %>%
    rename(nuclear = paste0(split, "_nuclear"),
           cytosolic = paste0(split, "_cytosolic"))

  lengths_list[[split]] <-
    lengths %>%
    select(
      transcript_id, whole_cell,
      paste0(split, "_nuclear"), paste0(split, "_cytosolic")
    ) %>%
    rename(nuclear = paste0(split, "_nuclear"),
           cytosolic = paste0(split, "_cytosolic"))

  n_molecules_list[[split]] <-
    n_molecules %>%
    select(
      transcript_id, whole_cell,
      paste0(split, "_nuclear"), paste0(split, "_cytosolic")
    ) %>%
    rename(nuclear = paste0(split, "_nuclear"),
           cytosolic = paste0(split, "_cytosolic"))
}
