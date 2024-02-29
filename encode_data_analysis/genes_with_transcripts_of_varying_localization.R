#===========================================================
# How many genes with transcripts of different localization?
#===========================================================

# per cell line
pdf("figures/encode_data/final/description_of_strongly_localized_transcripts/genes_with_transcripts_of_varying_localization.pdf",
  family = font,
  paper = "letter",
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
  mutate(gene_id = annotation$gene_id[idxs]) %>%
  mutate(transcript_type = annotation$transcript_type[idxs]) %>%
  mutate(`< 0.4` = ifelse(loc_idx < 0.4, TRUE, FALSE)) %>%
  mutate(`> 0.6` = ifelse(loc_idx > 0.6, TRUE, FALSE)) %>%
  mutate(middle = !(`< 0.4`) & !(`> 0.6`)) %>%
  group_by(cell_line, gene_id) %>%
  summarise(
    total = n(), cytosolic = sum(`> 0.6`),
    nuclear = sum(`< 0.4`), middle = sum(middle),
    n_trans_type =
      transcript_type %>% unique() %>% length()
  ) %>%
  mutate(max = pmax(cytosolic, nuclear, middle)) %>%
  summarise(
    `# genes expressed` = n(),
    `# transcripts expressed` = sum(total),
    `# cytosolic transcripts` = sum(cytosolic),
    `% cytosolic transcripts` = 100 * `# cytosolic transcripts` / `# transcripts expressed`,
    `# nuclear transcripts` = sum(nuclear),
    `% nuclear transcripts` = 100 * `# nuclear transcripts` / `# transcripts expressed`,
    `# genes w/ > 1 expressed transcript` =
      sum(total > 1),
    `# genes w/ transcripts of different localization` =
      sum(total != max),
    `% genes w/ transcripts of different localization` =
      100 * sum(total != max) /
        `# genes w/ > 1 expressed transcript`,
    `# genes w/ transcripts of the same biotype and different localization` =
      sum((total != max) & (n_trans_type == 1)),
    `% genes w/ transcripts of the same biotype and different localization` =
      100 * sum((total != max) & (n_trans_type == 1)) /
        `# genes w/ > 1 expressed transcript`,
  ) %>%
  mutate(across(starts_with("%"), ~ round(.x, 1))) %>%
  mutate(across(everything(), ~ format(.x, big.mark = ","))) %>%
  pivot_longer(-1,
    names_to = "measurement",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = cell_line, values_from = value) %>%
  mutate(measurement = str_wrap(measurement, 20)) %>%
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
