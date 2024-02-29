#============================================================================
# Barplot showing the proportion of transcripts that are cytosolic or nuclear
# in all cell lines where they are expressed
#============================================================================

pdf(
  "figures/encode_data/final/description_of_strongly_localized_transcripts/proportions_of_nuc_cyt_transcripts.pdf",
  width = 3.5,
  height = 3
)

loc_n_cell %>%
  filter(cutoff %in% c("<= 0.5", ">= 0.5") & n_cell == n_expressed) %>%
  select(-c(n_cell)) %>%
  group_by(cutoff, n_expressed) %>%
  summarise(pass = n()) %>%
  left_join(
    loc_n_cell %>%
      select(transcript_id, n_expressed) %>%
      distinct() %>%
      group_by(n_expressed) %>%
      summarise(total = n())
  ) %>%
  mutate(
    prop = pass / total,
    loc = ifelse(cutoff == "<= 0.5", "nuclear", "cytosolic")
  ) %>%
  ungroup() %>%
  select(-c(total, cutoff)) %>%
  mutate(
    n_expressed = factor(n_expressed),
    loc = factor(loc, levels = c("nuclear", "cytosolic"))
  ) %>%
  ggplot(aes(x = n_expressed, y = prop, fill = loc)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    fill = "Localization",
    x = "Number of cell lines",
    y = "Proportion of transcripts",
    tag = "C"
  ) +
  geom_text(aes(label = pass),
    position = position_stack(vjust = 0.5),
    size = 2.8, family = font
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text = element_text(size = size_small),
    legend.position = "bottom"
  )

dev.off()
