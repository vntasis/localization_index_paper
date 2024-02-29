#============================================
## Association of the localization index with
## different transcript structure features
#============================================

loc_idxs_final %<>%
  left_join(transcript_str)


# Use the transcripts with consistent localization in all cell lines
x <-
  loc_n_cell %>%
  ungroup() %>%
  filter(n_cell == n_expressed &
    cutoff %in% c("> 0.6", "< 0.4") & n_cell > 4) %>%
  mutate(loc = ifelse(cutoff == "< 0.4", "nuclear", "cytosolic")) %>%
  select(transcript_id, loc) %>%
  left_join(loc_idxs_final) %>%
  select(-c(loc_idx, cell_line)) %>%
  distinct() %>%
  select(
    loc, pct_gc_exons, total_exon_len,
    total_intron_len, nb_introns, "5_utr_len", "3_utr_len"
  ) %>%
  pivot_longer(2:7, "feature") %>%
  mutate(loc = factor(loc, levels = c("nuclear", "cytosolic"))) %>%
  mutate(
    feature =
      str_replace_all(
        feature,
        "nb_introns",
        "# Introns"
      ) %>%
      str_replace_all(
        "pct_gc_exons",
        "GC%"
      ) %>%
      str_replace_all(
        "total_exon_len",
        "(|Exon|)"
      ) %>%
      str_replace_all(
        "total_intron_len",
        "(|Intron|)"
      ) %>%
      str_replace_all(
        "5_utr_len",
        "|5' UTR|"
      ) %>%
      str_replace_all(
        "3_utr_len",
        "|3' UTR|"
      )
  )

meds <-
  x %>%
  group_by(loc, feature) %>%
  summarise(med = median(value, na.rm = TRUE)) %>%
  ungroup()

pdf("figures/encode_data/final/transcript_feature_association/localization_and_transcript_features_association.pdf",
  width = 4,
  height = 3.2
)

x %>%
  ggplot(aes(loc, value, colour = loc)) +
  geom_boxplot(fill = "grey") +
  facet_wrap(~feature, ncol = 3, scales = "free_y") +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text.y = element_text(size = size_small),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = size_big, face = "bold"),
    strip.background = element_rect(
      colour = "white",
      fill = "white"
    )
  ) +
  scale_y_log10(oob = scales::squish_infinite) +
  scale_colour_brewer(palette = "Dark2") +
  geom_text(
    data = meds, aes(y = med, label = round(med, 3)),
    size = 2.8, vjust = -0.5, family = font
  ) +
  labs(
    colour = "Localization",
    y = "", x = "",
    tag = "C"
  ) +
  stat_compare_means(
    label = "p.signif", method = "wilcox.test",
    comparisons = list(c("nuclear", "cytosolic")),
    vjust = 2, family = font, size = 2.8
  )


dev.off()
