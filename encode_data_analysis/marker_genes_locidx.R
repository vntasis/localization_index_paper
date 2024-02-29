#==============================================
## Localization index of all the transcripts of
## some example (marker) genes
#==============================================

markers <- c("GAPDH", "NEAT1", "NORAD", "CDK1")

pdf("figures/encode_data/final/locidx_in_transcript_examples/markers_SK-N-SH.pdf",
  width = 4,
  height = 2.5
)
loc_idxs_final %>%
  left_join(annotation_v45 %>%
    select(gene_id, gene_name) %>%
    distinct(), by = "gene_id") %>%
  filter(cell_line == "SK-N-SH" & gene_name %in% markers) %>%
  arrange(gene_name, transcript_type) %>%
  mutate(transcript_id = factor(transcript_id,
    levels = rev(transcript_id)
  )) %>%
  ggplot(aes(y = loc_idx, x = transcript_id, colour = gene_name)) +
  geom_point(aes(shape = transcript_type), size = 2) +
  geom_segment(aes(x = transcript_id, xend = transcript_id, y = 0,
                   yend = loc_idx), size = 1) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text.y = element_text(size = size_small),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    y = "Localization index (LI)",
    x = "Transcript",
    tag = "B",
    shape = "Transcript type"
  ) +
  scale_color_paletteer_d("ggsci::default_nejm", name = "Gene")
dev.off()
