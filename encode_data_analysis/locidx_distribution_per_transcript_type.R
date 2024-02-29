#************************************
## Distribution of localization index
## per transcript type
#************************************

loc_idxs_final %<>%
  separate(Transcript, c("transcript_id", "x"), sep = "[.]", remove = FALSE) %>%
  select(-x) %>%
  left_join(
    annotation_v45 %>%
      separate(transcript_id, c("transcript_id", "x"), sep = "[.]") %>%
      select(gene_id, gene_type, transcript_id, transcript_type)
  ) %>%
  mutate(transcript_id = Transcript) %>%
  select(-Transcript)

idxs <-
  loc_idxs_final %$%
  match(transcript_id, annotation$transcript_id)


# pooled
est_per_transcript_type <-
  loc_idxs_final %>%
  group_by(transcript_type) %>%
  summarise(
    number = n(), mean = mean(loc_idx),
    median = median(loc_idx), sd = sd(loc_idx)
  ) %>%
  arrange(desc(number))

meds <-
  loc_idxs_final %>%
  filter(transcript_type
  %in% c(
      "retained_intron", "nonsense_mediated_decay",
      "lncRNA", "protein_coding"
    )) %>%
  group_by(transcript_type) %>%
  summarise(med = median(loc_idx)) %>%
  ungroup() %>%
  mutate(
    transcript_type =
      str_replace_all(
        transcript_type,
        "retained_intron",
        "retained intron"
      )
  ) %>%
  mutate(
    transcript_type =
      str_replace_all(
        transcript_type,
        "protein_coding",
        "protein coding"
      )
  ) %>%
  mutate(
    transcript_type =
      str_replace_all(
        transcript_type,
        "nonsense_mediated_decay",
        "nonsense mediated decay"
      )
  )


pdf("figures/encode_data/final/locidx_distribution_per_transcript_type/pooled_violin_plots.pdf",
  width = 3.5,
  height = 2.5
)

loc_idxs_final %>%
  filter(transcript_type
  %in% c(
      "retained_intron", "nonsense_mediated_decay",
      "lncRNA", "protein_coding"
    )) %>%
  mutate(
    transcript_type =
      str_replace_all(
        transcript_type,
        "retained_intron",
        "retained intron"
      )
  ) %>%
  mutate(
    transcript_type =
      str_replace_all(
        transcript_type,
        "protein_coding",
        "protein coding"
      )
  ) %>%
  mutate(
    transcript_type =
      str_replace_all(
        transcript_type,
        "nonsense_mediated_decay",
        "nonsense mediated decay"
      )
  ) %>%
  ggplot(aes(transcript_type, loc_idx, colour = transcript_type)) +
  geom_violin(trim = TRUE, fill = "grey") +
  geom_boxplot(width = 0.05, fill = "grey") +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text.x = element_text(
      size = size_small,
      angle = 45, vjust = 1, hjust = 1
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_text(
    data = meds, aes(y = med, label = round(med, 3)),
    size = 2.8, vjust = -0.5, family = font
  ) +
  labs(
    x = "Transcript type",
    y = "Localization index",
    tag = "A"
  ) +
  coord_flip() +
  scale_colour_paletteer_d("ggsci::dark_uchicago",
    guide = guide_legend(reverse = TRUE)
  )
dev.off()


# per cell-line
meds <-
  loc_idxs_final %>%
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
  filter(transcript_type
  %in% c(
      "retained_intron", "nonsense_mediated_decay",
      "lncRNA", "protein_coding", "protein_coding_CDS_not_defined",
      "processed_pseudogene", "unprocessed_pseudogene",
      "transcribed_processed_pseudogene",
      "transcribed_unprocessed_pseudogene"
    )) %>%
  group_by(cell_line, transcript_type) %>%
  summarise(med = median(loc_idx)) %>%
  ungroup()


pdf(
  "figures/encode_data/final/locidx_distribution_per_transcript_type/violin_plots_per_cell_line.pdf",
  height = 8
)

loc_idxs_final %>%
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
  filter(transcript_type
  %in% c(
      "retained_intron", "nonsense_mediated_decay",
      "lncRNA", "protein_coding", "protein_coding_CDS_not_defined",
      "processed_pseudogene", "unprocessed_pseudogene",
      "transcribed_processed_pseudogene",
      "transcribed_unprocessed_pseudogene"
    )) %>%
  mutate(transcript_type = factor(transcript_type,
    levels =
      c(
        "lncRNA", "nonsense_mediated_decay",
        "protein_coding", "retained_intron",
        "protein_coding_CDS_not_defined",
        "processed_pseudogene",
        "unprocessed_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unprocessed_pseudogene"
      )
  )) %>%
  ggplot(aes(transcript_type, loc_idx, colour = transcript_type)) +
  geom_violin(trim = TRUE, fill = "grey") +
  geom_boxplot(width = 0.05, fill = "grey") +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    axis.text = element_text(size = size_small),
    legend.position = "none",
    strip.text = element_text(size = size_big, face = "bold"),
    strip.background = element_rect(
      colour = "white",
      fill = "white"
    )
  ) +
  geom_text(
    data = meds, aes(y = med, label = round(med, 3)),
    size = 2.8, vjust = -0.5, family = font
  ) +
  xlab("Transcript type") +
  ylab("Localization index") +
  facet_wrap(~cell_line, ncol = 3) +
  coord_flip() +
  scale_colour_paletteer_d("ggsci::dark_uchicago")
dev.off()


# comparison retained intron transcripts from pc and lnc genes
loc_idxs_final %>%
  filter(transcript_type == "retained_intron") %>%
  group_by(gene_type, cell_line) %>%
  summarise(
    number = n(), mean = mean(loc_idx),
    median = median(loc_idx), sd = sd(loc_idx)
  ) %>%
  arrange(cell_line, desc(number))

meds <-
  loc_idxs_final %>%
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
  filter(transcript_type == "retained_intron") %>%
  group_by(gene_type, cell_line) %>%
  summarise(med = median(loc_idx)) %>%
  ungroup()

pdf(
  "figures/encode_data/final/locidx_distribution_per_transcript_type/retained_intron_transcripts.pdf",
  height = 5,
  width = 5
)

loc_idxs_final %>%
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
  filter(transcript_type == "retained_intron") %>%
  ggplot(aes(gene_type, loc_idx, colour = gene_type)) +
  geom_violin(trim = TRUE, fill = "grey") +
  geom_boxplot(width = 0.05, fill = "grey") +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    axis.text = element_text(size = size_small),
    legend.position = "none",
    strip.text = element_text(size = size_big, face = "bold"),
    strip.background = element_rect(
      colour = "white",
      fill = "white"
    )
  ) +
  geom_text(
    data = meds, aes(y = med, label = round(med, 3)),
    size = 2.8, vjust = -0.5, family = font
  ) +
  xlab("Gene type") +
  ylab("Localization index") +
  facet_wrap(~cell_line, ncol = 3) +
  coord_flip() +
  scale_colour_paletteer_d("ggsci::dark_uchicago")
dev.off()
