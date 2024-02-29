#******************************************************************
## Turn the count (and lenths) tibble into a list of tibble.
## Each element of the list will correspond to
## whole cell, cytosolic, nuclear counts of a sample.
#******************************************************************

cell_types <- meta %>%
  pull("biosample") %>%
  unique()

#----------------------------
## Create the list of tibbles
#----------------------------


counts_list <- list()
lengths_list <- list()
np_idr_score_filter <- list()
for (cell_type in cell_types) {
  replicates <-
    meta %>%
    filter(biosample == cell_type) %$%
    replicate %>%
    unique()

  for (rep in replicates) {
    n_samples <-
      meta %>%
      filter(biosample == cell_type & replicate == rep) %>%
      nrow()

    if (n_samples < 3) next

    sample <- paste(cell_type, rep, sep = "_")

    counts_list[[sample]] <-
      counts %>%
      select(
        meta %>%
          filter(biosample == cell_type & replicate == rep) %$%
          file
      )

    lengths_list[[sample]] <-
      lengths %>%
      select(
        meta %>%
          filter(biosample == cell_type & replicate == rep) %$%
          file
      )

    idr_filter <- rep(TRUE, nrow(counts))

    for (i in 1:3) {
      id <- colnames(counts_list[[sample]])[i]
      fr <- meta %>%
        filter(file == id) %$%
        fraction
      colnames(counts_list[[sample]])[i] <- fr

      # npIDR related filter
      idr_filter %<>% `&`(np_idr_score[[id]] < 0.1)


      id <- colnames(lengths_list[[sample]])[i]
      fr <- meta %>%
        filter(file == id) %$%
        fraction
      colnames(lengths_list[[sample]])[i] <- fr
    }

    np_idr_score_filter[[sample]] <- idr_filter

    counts_list[[sample]] %<>%
      add_column(.before = 1, transcript_id = counts$transcript_id)
    lengths_list[[sample]] %<>%
      add_column(.before = 1, transcript_id = lengths$transcript_id)
  }
}

# Add A549, keratinocytes, MCF-7, endothelial
combine_reps <- function(cell_type, reps, n) {
  sample <- paste(cell_type, n, sep = "_")

  ids <- character()
  for (rep in reps) {
    ids <- c(ids,
             meta %>%
               filter(biosample == cell_type & replicate == rep) %$% file)
  }

  counts_list[[sample]] <<-
    counts %>%
    select(ids)

  lengths_list[[sample]] <<-
    lengths %>%
    select(ids)

  idr_filter <- rep(TRUE, nrow(counts))

  for (i in 1:3) {
    id <- colnames(counts_list[[sample]])[i]
    fr <- meta %>%
      filter(file == id) %$%
      fraction
    colnames(counts_list[[sample]])[i] <<- fr

    # npIDR related filter
    idr_filter %<>% `&`(np_idr_score[[id]] < 0.1)

    id <- colnames(lengths_list[[sample]])[i]
    fr <- meta %>%
      filter(file == id) %$%
      fraction
    colnames(lengths_list[[sample]])[i] <<- fr
  }

  np_idr_score_filter[[sample]] <<- idr_filter

  counts_list[[sample]] <<-
    counts_list[[sample]] %>%
    add_column(.before = 1, transcript_id = counts$transcript_id)
  lengths_list[[sample]] <<-
    lengths_list[[sample]] %>%
    add_column(.before = 1, transcript_id = lengths$transcript_id)
}
combine_reps("A549", c(1, 3), "1")
combine_reps("A549", c(2, 4), "2")
combine_reps("keratinocyte", c(1, 3), "1")
combine_reps("keratinocyte", c(2, 4), "2")
combine_reps("MCF-7", c(1, 3), "1")
combine_reps("MCF-7", c(2, 4), "2")
combine_reps("endothelial", c(1, 3), "1")
combine_reps("endothelial", c(2, 4), "2")

rm(list = c("cell_type", "replicates", "rep", "sample",
            "id", "fr", "i", "n_samples", "combine_reps", "idr_filter"))
