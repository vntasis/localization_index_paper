#-------------------------
## Load required libraries
#-------------------------

library(tidyverse)
library(magrittr)
library(vroom)
library(edgeR)
library(grid)
library(gridExtra)
library(gtable)

#---------------------------------------------
## Load Data:
## - raw counts (expected counts from RSEM)
## - effective lengths
## - Annotation data
## - Metadata of the single-cell data
#---------------------------------------------

counts <- vroom("../SINCseq_data/expected_counts_CNT.tsv")
lengths <- vroom("../SINCseq_data/lengths_CNT.tsv")
lengths %<>% slice(match(counts$transcript_id, lengths$transcript_id))

annotation <- vroom::vroom("../Data/annotation_transcripts.tsv",
  "\t",
  col_names = c(
    "chr", "start", "end",
    "strand", "gene_id",
    "transcript_id",
    "transcript_type",
    "gene_name"
  )
)

meta <- vroom("../SINCseq_data/meta.tsv")

#-----------------------------------------
## Create a table with the metadata in pdf
#-----------------------------------------

size_small <- 8
size_medium <- 9
size_big <- 10
frame_size <- 2
font <- "NimbusSan"


pdf("figures/sincseq/data_accession.pdf",
  family = font,
  paper = "letter",
  height = 9
)
meta %>%
  select(Run, Experiment, Sample, LibraryName) %>%
  arrange(LibraryName) %>%
  tableGrob(
    rows = NULL,
    theme =
      ttheme_default(
        base_size = size_small,
        padding = unit(c(3, 3), "mm")
      )
  ) %>%
  gtable_add_grob(
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 2, b = nrow(.), l = 1, r = ncol(.)
  ) %>%
  gtable_add_grob(
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 1, l = 1, r = ncol(.)
  ) %>%
  grid.draw()
dev.off()

#-----------------
## Preprocess data
#-----------------

files <- c("scripts/sinc_seq_data/preprocess.R")

sapply(files, source, echo = TRUE)
