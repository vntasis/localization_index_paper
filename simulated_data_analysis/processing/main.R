#-------------------------
## Load required libraries
#-------------------------

library(tidyverse)
library(magrittr)
library(vroom)
library(edgeR)
library(paletteer)
library(ggpointdensity)


#---------------------------------------------
## Load Data:
## - raw counts (expected counts from RSEM)
## - effective lengths
## - simulated number of molecules per sample
#---------------------------------------------

counts <- vroom("TranscriptQuantifications_raw.tsv")
lengths <- vroom("TranscriptQuantifications_ef_length.tsv")
n_molecules <- vroom("Nmolecules.tsv")


splits <- paste0("split", seq(50, 80, 10))

#-----------------
## Preprocess data
#-----------------

files <- c(
  "scripts/simulation_analysis/restructure_counts.R",
  "scripts/simulation_analysis/preprocess.R"
)

sapply(files, source, echo = TRUE)


#-------------------------------
## Save fpkm values in a R image
#-------------------------------
save(
  list = c("fpkms", "counts_list", "lengths_list", "n_molecules_list"),
  file = "RData/simulated_data_final.RData"
)



#*************************************************
## RUN STAN-NF with our model for esitmating beta*
#*************************************************


#--------------------------
## Load estimations of beta
#--------------------------
load("RData/posterior_estimations_simulation_data.RData")



#-------------------------
## Set variables for plots
#-------------------------
size_small <- 8
size_medium <- 9
size_big <- 10
frame_size <- 2
font <- "NimbusSan"

#------------------------------------------------------------------------
## Calculate for every split the loclalizaiton index
## for all transcripts based on:
## - The actual (cytosolic / whole_cell) number of molecules (actual)
## - The (cytosolic / (nuclear + cytosolic)) fpkms (estimated)
## - The same ratio corrected by our beta estimates (corrected)
#------------------------------------------------------------------------

loc_idx <- tibble(
  transcript_id = character(),
  split = character(),
  method = character(),
  loc_idx = numeric()
)

for (split in splits) {
  loc_idx %<>%
    bind_rows(
              n_molecules_list[[split]] %>%
                mutate(loc_idx = cytosolic / whole_cell) %>%
                drop_na() %>%
                select(transcript_id, loc_idx) %>%
                add_column(split = split, .before = 2) %>%
                add_column(method = "actual", .before = 3)) %>%
    bind_rows(
              fpkms[[split]] %>%
                mutate(loc_idx = (cytosolic / (cytosolic + nuclear))) %>%
                select(transcript_id, loc_idx) %>%
                add_column(split = split, .before = 2) %>%
                add_column(method = "estimated", .before = 3)) %>%
    bind_rows(
              fpkms[[split]] %>%
                add_column(beta = posterior[[split]][[1]]) %>%
                mutate(loc_idx = (beta * cytosolic
                                  / (beta * cytosolic
                                     + (1 - beta) * nuclear))) %>%
                select(transcript_id, loc_idx) %>%
                add_column(split = split, .before = 2) %>%
                add_column(method = "corrected", .before = 3))
}

loc_idx %<>%
  mutate(
    split =
      replace(split, split == "split50", "0.5") %>%
      replace(., . == "split60", "0.6") %>%
      replace(., . == "split70", "0.7") %>%
      replace(., . == "split80", "0.8")
  )

## Plot loc_idx distribtution
pdf("figures/simulated_data/loc_idx_simulated_data.pdf",
  height = 4,
  width = 4
)
loc_idx %>%
  mutate(method = factor(method,
    levels = c("actual", "estimated", "corrected"),
    labels = c("Simulated C/W", "Naive LI", "LI")
  )) %>%
  ggplot(aes(x = loc_idx)) +
  geom_histogram(fill = "#a6bddb", col = "black", size = 0.2) +
  facet_grid(vars(method), vars(split), scales = "free_y") +
  scale_x_continuous(n.breaks = 3) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    strip.text = element_text(size = size_big, face = "bold"),
    strip.background =
      element_rect(colour = "white", fill = "white"),
    plot.tag = element_text(face = "bold"),
    axis.text = element_text(size = size_small)
  ) +
  labs(
    y = "Frequency",
    x = "",
    tag = "A"
  )
dev.off()



#--------------------------------------------------
## For every simulated transcript calculated
## the error of the estimation of loc_idx
## both for corrected and not corrected estimations
#--------------------------------------------------

pdf("figures/simulated_data/error_loc_idx_simulated_data.pdf",
  width = 2, height = 3
)
loc_idx %>%
  pivot_wider(names_from = method, values_from = loc_idx) %>%
  drop_na() %>%
  mutate(
    `Naive LI` = estimated - actual,
    LI = corrected - actual
  ) %>%
  select(transcript_id, split, `Naive LI`, LI) %>%
  pivot_longer(!c("transcript_id", "split"),
    values_to = "error",
    names_to = "error_type"
  ) %>%
  ggplot(aes(split, error, colour = error_type)) +
  geom_boxplot(outlier.size = .1, lwd = .1, fill = "grey") +
  facet_wrap(~error_type, ncol = 1) +
  scale_color_manual(values = rev(paletteer_d("ggsci::dark_uchicago")[1:2])) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    strip.text = element_text(size = size_big, face = "bold"),
    strip.background =
      element_rect(colour = "white", fill = "white"),
    legend.position = "none",
    plot.tag = element_text(face = "bold"),
    axis.text = element_text(size = size_small)
  ) +
  labs(
    y = "Error",
    x = "Simulated beta",
    tag = "B"
  )
dev.off()


## Check the correlation of the esimation error with the actual loc_idx
loc_idx_cor <-
  loc_idx %>%
  pivot_wider(names_from = method, values_from = loc_idx) %>%
  drop_na() %>%
  mutate(
    Not_corrected = estimated - actual,
    Corrected = corrected - actual
  ) %>%
  select(transcript_id, split, actual, Not_corrected, Corrected) %>%
  group_by(split) %>%
  summarise(
    Not_corrected = cor(actual, Not_corrected, method = "spearman"),
    Corrected = cor(actual, Corrected, method = "spearman")
  ) %>%
  pivot_longer(!split, names_to = "error_type", values_to = "cor") %>%
  mutate(cor = paste0("r = ", round(cor, 3)))

pdf("figures/simulated_data/cor_idx_error.pdf",
  height = 8, width = 6
)
loc_idx %>%
  pivot_wider(names_from = method, values_from = loc_idx) %>%
  drop_na() %>%
  mutate(
    Not_corrected = estimated - actual,
    Corrected = corrected - actual
  ) %>%
  select(transcript_id, split, actual, Not_corrected, Corrected) %>%
  pivot_longer(!c("transcript_id", "split", "actual"),
    values_to = "error",
    names_to = "error_type"
  ) %>%
  ggplot(aes(actual, error)) +
  geom_pointdensity() +
  geom_smooth(method = MASS::rlm, color = "red") +
  facet_grid(vars(split), vars(error_type)) +
  geom_text(
    data = loc_idx_cor,
    aes(x = 0.3, y = 1, label = cor),
    size = 4, vjust = 0.5
  ) +
  scale_color_paletteer_c("viridis::viridis") +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    axis.text = element_text(size = size_small)
  )
dev.off()
