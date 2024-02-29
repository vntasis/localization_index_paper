#-------------------------
## Load required libraries
#-------------------------

library(tidyverse)
library(magrittr)
library(vroom)
library(edgeR)
library(paletteer)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)

options(scipen = 999)

#---------------------------------------------
## Load Data:
## - raw counts (expected counts from RSEM)
## - effective lengths
## - Annotation data
## - Transcript str info (GC%, length, etc.)
#---------------------------------------------

counts <- vroom("../encode_data/transcript_data/counts.tsv")
lengths <- vroom("../encode_data/transcript_data/lengths.tsv")
meta <- vroom("../encode_data/transcript_data/meta.tsv")
meta %<>%
  rename(
    file = "File accession",
    experiment = "Experiment accession",
    biosample = "Biosample term name",
    replicate = "Biological replicate(s)",
    fraction = "Fraction"
  )

meta %<>%
  mutate(biosample
         = replace(biosample,
                   biosample == "endothelial cell of umbilical vein",
                   "endothelial"))

annotation <-
  vroom::vroom("../encode_data/transcript_data/transcript_annotation.tsv",
    "\t",
    col_names = c(
      "chr", "start", "end",
      "strand", "gene_id",
      "gene_type",
      "transcript_id",
      "transcript_type",
      "gene_name"
    )
  )

annotation_v45 <-
  vroom::vroom("../encode_data/transcript_data/transcript_annotation_v45.tsv",
    "\t",
    col_names = c(
      "chr", "start", "end",
      "strand", "gene_id",
      "gene_type",
      "transcript_id",
      "transcript_type",
      "gene_name"
    )
  )

transcript_str <- vroom("../encode_data/transcript_data/transcript_str.tsv")

#-----------------------------------------------------
## Preprocess data and load necessary custom functions
#-----------------------------------------------------

files <- c(
  "scripts/encode_data/apply_npidr.R",
  "scripts/encode_data/restructure_counts.R",
  "scripts/encode_data/preprocess.R"
)

sapply(files, source, echo = TRUE)


#----------------------------
## Estimate beta with Stan-NF
#----------------------------
# Load estimates
load("RData/posterior_estimations_encode_data_transcripts_npIDR.RData")


#-----------------------------------
## Summarise the different estimates
#-----------------------------------


estimated <- tibble(sample = character(), beta = numeric())
for (sample in names(posterior)) {
  estimated %<>% add_row(sample = sample, beta = posterior[[sample]][["beta"]])
}

estimated %<>%
  add_column(biosample
             = str_split(.$sample, "_", n = 2, simplify = TRUE) %>%
               `[`(, 1))

estimated %<>%
  group_by(biosample) %>%
  summarise(beta_med = median(beta)) %>%
  inner_join(estimated) %>%
  arrange(desc(beta_med), desc(beta)) %>%
  select(-beta_med) %>%
  mutate(sample = factor(sample, levels = sample))

size_small <- 8
size_medium <- 9
size_big <- 10
frame_size <- 2
font <- "NimbusSan"

pdf("figures/encode_data/final/beta_and_locidx_estimates/estimated_betas.pdf",
  width = 3.2,
  height = 3
)
estimated %>%
  mutate(
    biosample =
      str_replace_all(
        biosample,
        "endothelial",
        "HUVEC"
      )
  ) %>%
  mutate(
    biosample =
      str_replace_all(
        biosample,
        "keratinocyte",
        "NHEK"
      )
  ) %>%
  ggplot(aes(x = sample, y = beta, fill = biosample)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text.y = element_text(size = size_small),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.box.margin = margin(t = -1, unit = "in"),
    plot.margin = margin(t = 1, unit = "in")
  ) +
  labs(
    y = "Estimated beta",
    x = "",
    tag = "A"
  ) +
  scale_fill_paletteer_d("ggsci::default_igv", name = "ENCODE sample")
dev.off()

samples_selected <-
  estimated %>%
  filter(beta <= 0.95) %>%
  pull(sample) %>%
  as.character()


#--------------------------------------------
## Check localization index across cell lines
#--------------------------------------------

# Function that calcules the localization indices for a sample
get_loc_idx <- function(samp) {
  beta <-
    estimated %>%
    filter(sample == samp) %>%
    pull(beta)

  fpkms[[samp]] %>%
    add_column(
      cyt_prop = (.$cytosolic * beta)
      / (.$cytosolic * beta + .$nuclear * (1 - beta))
    )
}

# Function that calcules the naive localization indices for a sample
get_naive_loc_idx <- function(samp) {
  beta <-
    estimated %>%
    filter(sample == samp) %>%
    pull(beta)

  fpkms[[samp]] %>%
    add_column(
      cyt_prop = (.$cytosolic) / (.$cytosolic + .$nuclear)
    )
}

# Run the above function for every sample in the dataset
loc_idxs <- list()
for (s in samples_selected) {
  loc_idxs[[s]] <-
    get_loc_idx(s) %>%
    select(transcript_id, cyt_prop) %>%
    rename(!!s := cyt_prop)
}

# Merge everything in one table
loc_idxs_tbl <- reduce(loc_idxs, full_join)

loc_idxs_tbl %<>%
  pivot_longer(-1, names_to = "Sample", values_to = "loc_idx") %>%
  rename(Transcript = transcript_id) %>%
  replace_na(list(loc_idx = -1))

loc_idxs_tbl %>%
  filter(loc_idx != -1) %>%
  group_by(Sample) %>%
  summarise(number = n())


# Average values across biological replicates
loc_idxs_tbl %<>%
  separate(Sample, c("cell_line", "x"),
           sep = "_", remove = FALSE) %>%
  select(-x)

cell_lines <-
  loc_idxs_tbl %>%
  pull(cell_line) %>%
  unique()

samps1rep <-
  c("GM12878", "H1", "HeLa-S3", "IMR-90", "K562")


loc_idxs_final <-
  loc_idxs_tbl %>%
  filter(cell_line %in% samps1rep) %>%
  filter(loc_idx != -1) %>%
  select(-Sample) %>%
  arrange(cell_line)

for (s in cell_lines) {
  if (s %in% samps1rep) next

  selected_transcripts <-
    loc_idxs_tbl %>%
    filter(cell_line == s) %>%
    filter(loc_idx != -1) %>%
    group_by(Transcript) %>%
    summarise(n = n()) %>%
    filter(n == 2) %>%
    pull(Transcript)

  loc_idxs_final <-
    loc_idxs_tbl %>%
    filter(cell_line == s) %>%
    filter(loc_idx != -1) %>%
    filter(Transcript %in% selected_transcripts) %>%
    group_by(Transcript, cell_line) %>%
    summarise(loc_idx = mean(loc_idx)) %>%
    ungroup() %>%
    bind_rows(loc_idxs_final)
}


# Make a heatmap ilustraging loc_idxs_tbl
library(tidyHeatmap)
pdf("figures/encode_data/final/beta_and_locidx_estimates/localization_index_heatmap.pdf",
  family = font
)
loc_idxs_final %>%
  filter(loc_idx != -1) %>%
  mutate(
    cell_line =
      str_replace_all(
        cell_line,
        "endothelial",
        "HUVEC"
      )
  ) %>%
  mutate(
    cell_line =
      str_replace_all(
        cell_line,
        "keratinocyte",
        "NHEK"
      )
  ) %>%
  pivot_wider(names_from = cell_line, values_from = loc_idx) %>%
  drop_na() %>%
  pivot_longer(-1, names_to = "Cell line", values_to = "loc_idx") %>%
  heatmap(Transcript, `Cell line`, loc_idx,
    .scale = "none",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    show_column_names = FALSE,
    palette_value =
      circlize::colorRamp2(
        seq(0, 1, length.out = 11),
        paletteer_d("RColorBrewer::RdBu")
      )
  ) %>%
  add_tile(`Cell line`, paletteer_d("ggsci::default_igv")[2:11])
dev.off()

# Clean up
rm(list = c("s", "selected_transcripts", "sample"))



#-------------------------------------------------
## Study the association of the localization index
## with different transcript features
#-------------------------------------------------
source("scripts/encode_data/feature_association.R", echo = TRUE)
