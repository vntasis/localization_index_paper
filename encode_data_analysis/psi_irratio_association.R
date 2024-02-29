#============================================
## Association of exon PSI and intron IRratio
## with the transcript localization
#============================================

# Load data and metadata
meta_bam <- vroom("../encode_data/bam_files/meta.tsv",
  col_names = c(
    "file", "experiment", "biosample",
    "replicate", "fraction"
  ),
  skip = 1
)

meta_bam %<>%
  mutate(
    biosample =
      str_replace_all(
        biosample,
        "endothelial cell of umbilical vein",
        "endothelial"
      )
  )

exon_psi <- vroom("../encode_data/splicing_analysis/exon_psi.tsv")


# Make table with the accession numbers of the samples used
# For supplementary data
pdf("figures/encode_data/final/data_accession.pdf",
  family = font,
  paper = "letter",
  height = 10
)
meta %>%
  left_join(meta_bam,
    by = c(
      "experiment", "biosample",
      "replicate", "fraction"
    )
  ) %>%
  rename(quantification_file = file.x, bam_file = file.y) %>%
  select(1, 6, 2, 3, 4, 5) %>%
  filter(!(biosample == "H1" & replicate == "1")) %>%
  filter(!(biosample == "keratinocyte" & replicate == "5")) %>%
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
  tableGrob(
    rows = NULL,
    theme =
      ttheme_default(
        base_size = size_small,
        padding = unit(c(1.5, 1.5), "mm")
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


# Preprocess PSI data
exon_psi %<>%
  unite("exon", chr:strand, sep = "|") %>%
  pivot_longer(-c(1, 2), names_to = "file", values_to = "psi") %>%
  left_join(meta_bam) %>%
  filter(biosample %in% cell_lines) %>%
  rename(cell_line = biosample) %>%
  unite("biosample", cell_line:replicate, sep = "_", remove = FALSE) %>%
  filter(biosample %in% samples_selected) %>%
  select(-biosample)

# Average PSI per exon across the two replicates
exon_psi_final <-
  exon_psi %>%
  filter(cell_line %in% samps1rep) %>%
  select(-c("fraction", "file", "experiment", "replicate")) %>%
  arrange(cell_line)

for (s in cell_lines) {
  if (s %in% samps1rep) next

  exon_psi_final <-
    exon_psi %>%
    filter(cell_line == s) %>%
    select(-c("fraction", "file", "experiment", "replicate")) %>%
    group_by(exon, transcript_id, cell_line) %>%
    summarise(psi = mean(psi, na.rm = TRUE)) %>%
    ungroup() %>%
    bind_rows(exon_psi_final)
}

# Merge PSI data with the localization index info
exon_psi_final <-
  loc_idxs_final %>%
  select("transcript_id", "cell_line", "loc_idx") %>%
  left_join(exon_psi_final) %>%
  drop_na()

# Get nuclear and cytosolic frequency for every exon
exon_frequency <-
  exon_psi_final %>%
  select(-c("transcript_id", "psi")) %>%
  add_column(loc_idx_group = "mid") %>%
  mutate(loc_idx_group = replace(loc_idx_group, loc_idx > 0.6, "high")) %>%
  mutate(loc_idx_group = replace(loc_idx_group, loc_idx < 0.4, "low")) %>%
  group_by(exon, loc_idx_group, cell_line) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = loc_idx_group, values_from = count) %>%
  mutate(
    low = replace_na(low, 0),
    mid = replace_na(mid, 0),
    high = replace_na(high, 0)
  )


exon_frequency %>%
  filter(low == 0 & mid == 0) %$%
  table(high)



## Calculate the enrichment of
## low psi exons in nuclear and
## cytosolic specific exons
## ============================

# Cytosolic specific exons
cytosolic_exons <-
  exon_frequency %>%
  filter(low == 0 & mid == 0) %>%
  select(exon, cell_line)

N_cytosolic_exons <-
  cytosolic_exons %>%
  group_by(cell_line) %>%
  summarise(tot_cyt = n())

# Nuclear specific exons
nuclear_exons <-
  exon_frequency %>%
  filter(high == 0 & mid == 0) %>%
  select(exon, cell_line)

N_nuclear_exons <-
  nuclear_exons %>%
  group_by(cell_line) %>%
  summarise(tot_nuc = n())



n_low_psi <-
  exon_psi_final %>%
  group_by(cell_line) %>%
  group_map(function(x, y) {
    c <- y %>% pull(cell_line)
    x %>%
      filter(exon %in%
        (nuclear_exons %>%
          filter(cell_line == c) %>%
          pull(exon))) %>%
      select(-c("loc_idx", "transcript_id")) %>%
      distinct() %>%
      summarise(nuc = sum(psi < 0.5)) %>%
      add_column(
        cyt =
          (x %>%
            filter(exon %in%
              (cytosolic_exons %>%
                filter(cell_line == c) %>%
                pull(exon))) %>%
            select(-c("loc_idx", "transcript_id")) %>%
            distinct() %$%
            sum(psi < 0.5))
      ) %>%
      add_column(cell_line = c)
  }) %>%
  reduce(bind_rows)


# Number of low and high psi exons included in
# the whole set of exons in our dataset
N_low_psi <-
  exon_psi_final %>%
  select(-c("loc_idx", "transcript_id")) %>%
  distinct() %>%
  group_by(cell_line) %>%
  summarise(
    low = sum(psi < 0.5),
    high = sum(psi >= 0.5),
    total = n()
  )


# Calculate an enrichment
# and a p-value for the enrichment
low_psi_enr <-
  n_low_psi %>%
  left_join(N_cytosolic_exons) %>%
  left_join(N_nuclear_exons) %>%
  left_join(N_low_psi) %>%
  mutate(
    cytosolic = (cyt / tot_cyt) / (low / total),
    nuclear = (nuc / tot_nuc) / (low / total)
  ) %>%
  group_by(cell_line) %>%
  summarise(
    nuc_pvalue =
      phyper((nuc - 1), low,
        high, tot_nuc,
        lower.tail = FALSE,
        log.p = FALSE
      ),
    cyt_pvalue =
      phyper((cyt - 1), low,
        high, tot_cyt,
        lower.tail = FALSE,
        log.p = FALSE
      ),
    nuclear = nuclear,
    cytosolic = cytosolic
  )

low_psi_enr %<>%
  mutate(
    cell_line =
      str_replace_all(
        cell_line,
        "endothelial",
        "HUVEC"
      ) %>%
        str_replace_all(
          "keratinocyte",
          "NHEK"
        )
  ) %>%
  pivot_longer(4:5, names_to = "frac", values_to = "enr") %>%
  pivot_longer(2:3, names_to = "frac_pval", values_to = "pvalue") %>%
  filter((frac == "nuclear" & frac_pval == "nuc_pvalue") |
    (frac == "cytosolic" & frac_pval == "cyt_pvalue")) %>%
  select(-frac_pval) %>%
  mutate(pvalue = -log10(p.adjust(pvalue, method = "bonferroni"))) %>%
  add_column(sig_level = "") %>%
  mutate(sig_level = replace(sig_level, pvalue >= -log10(0.05), "*"))


# Plot enrichment
# ================
pdf("figures/encode_data/final/exon_skipping_intron_inclusion_enrichment/low_psi_enrichment.pdf",
  width = 3.2,
  height = 3
)

low_psi_enr %>%
  ggplot(aes(y = enr, x = frac, fill = cell_line)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(
    aes(
      y = enr + 0.05, x = frac,
      label = sig_level, group = cell_line
    ),
    position = position_dodge(width = .9),
    size = 2.8, family = font
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text = element_text(size = size_small),
    axis.text.x =
      element_text(angle = 45, vjust = 1, hjust = 1),
    legend.box.margin = margin(t = -1, b = -0.5, unit = "in"),
    plot.margin = margin(t = 1, unit = "in")
  ) +
  labs(
    y = "Enrichment of low PSI exons",
    x = "Fraction specific exons",
    tag = "D"
  ) +
  scale_fill_manual(
    values = paletteer_d("ggsci::default_igv")[2:11],
    name = "ENCODE sample"
  )
dev.off()



# IRratio
# ========
ir_ratio <- vroom("../encode_data/transcript_data/IRratio.tsv")
intron_coords <-
  vroom("../encode_data/transcript_data/introns.v29.primary_assembly.bed",
    col_names = c("chr", "start", "end", "intron", "intron_num", "strand")
  )

intron_coords %<>%
  separate(intron, c("x", "gene_id", "transcript_id"), sep = "[|]") %>%
  unite("intron", c("transcript_id", "intron_num"), sep = "-") %>%
  unite("coords", c("chr", "start", "end", "strand"), sep = "|") %>%
  select(-c("x", "gene_id"))

# Preprocess IRratio data
ir_ratio %<>%
  unite("intron", transcript_id:intron_num,
    sep = "-", remove = FALSE
  ) %>%
  pivot_longer(-c(1:3),
    names_to = "file",
    values_to = "irratio"
  ) %>%
  left_join(meta_bam) %>%
  filter(biosample %in% cell_lines) %>%
  rename(cell_line = biosample) %>%
  unite("biosample", cell_line:replicate, sep = "_", remove = FALSE) %>%
  filter(biosample %in% samples_selected) %>%
  select(-biosample)

# Average IRratio per intron across the two replicates
ir_ratio_final <-
  ir_ratio %>%
  filter(cell_line %in% samps1rep) %>%
  select(-c("fraction", "file", "experiment", "replicate")) %>%
  arrange(cell_line)

for (s in cell_lines) {
  if (s %in% samps1rep) next

  ir_ratio_final <-
    ir_ratio %>%
    filter(cell_line == s) %>%
    select(-c("fraction", "file", "experiment", "replicate")) %>%
    group_by(intron, transcript_id, intron_num, cell_line) %>%
    summarise(irratio = mean(irratio, na.rm = TRUE)) %>%
    ungroup() %>%
    bind_rows(ir_ratio_final)
}

# Merge IRratio data with the localization index info
ir_ratio_final <-
  loc_idxs_final %>%
  select("transcript_id", "cell_line", "loc_idx") %>%
  left_join(ir_ratio_final) %>%
  drop_na()

# Merge with intron coordinates
ir_ratio_final %<>%
  left_join(intron_coords) %>%
  select(-intron) %>%
  rename(intron = coords)

# Get nuclear and cytosolic frequency for every exon
intron_frequency <-
  ir_ratio_final %>%
  select(-c("transcript_id", "irratio", "intron_num")) %>%
  add_column(loc_idx_group = "mid") %>%
  mutate(loc_idx_group = replace(loc_idx_group, loc_idx > 0.6, "high")) %>%
  mutate(loc_idx_group = replace(loc_idx_group, loc_idx < 0.4, "low")) %>%
  group_by(intron, loc_idx_group, cell_line) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = loc_idx_group, values_from = count) %>%
  mutate(
    low = replace_na(low, 0),
    mid = replace_na(mid, 0),
    high = replace_na(high, 0)
  )



## Calculate the enrichment of
## high IRratio introns in nuclear
## and cytosolic specific introns
## ================================

# Cytosolic specific introns
cytosolic_introns <-
  intron_frequency %>%
  filter(low == 0 & mid == 0) %>%
  select(intron, cell_line)

N_cytosolic_introns <-
  cytosolic_introns %>%
  group_by(cell_line) %>%
  summarise(tot_cyt = n())

# Nuclear specific introns
nuclear_introns <-
  intron_frequency %>%
  filter(high == 0 & mid == 0) %>%
  select(intron, cell_line)

N_nuclear_introns <-
  nuclear_introns %>%
  group_by(cell_line) %>%
  summarise(tot_nuc = n())


# Number of high irratio introns included in
# nuclear and cytosolic specific introns

n_high_irratio <-
  ir_ratio_final %>%
  group_by(cell_line) %>%
  group_map(function(x, y) {
    c <- y %>% pull(cell_line)
    x %>%
      filter(intron %in%
        (nuclear_introns %>%
          filter(cell_line == c) %>%
          pull(intron))) %>%
      select(-c("loc_idx", "transcript_id", "intron_num")) %>%
      distinct() %>%
      summarise(nuc = sum(irratio > 0.5)) %>%
      add_column(
        cyt =
          (x %>%
            filter(intron %in%
              (cytosolic_introns %>%
                filter(cell_line == c) %>%
                pull(intron))) %>%
            select(-c("loc_idx", "transcript_id", "intron_num")) %>%
            distinct() %$%
            sum(irratio > 0.5))
      ) %>%
      add_column(cell_line = c)
  }) %>%
  reduce(bind_rows)


# Number of low and high irratio introns included in
# the whole set of introns in our dataset
N_irratio <-
  ir_ratio_final %>%
  select(-c("loc_idx", "transcript_id", "intron_num")) %>%
  distinct() %>%
  group_by(cell_line) %>%
  summarise(
    low = sum(irratio <= 0.5),
    high = sum(irratio > 0.5),
    total = n()
  )


# Calculate an enrichment
# and a p-value for the enrichment
high_irratio_enr <-
  n_high_irratio %>%
  left_join(N_cytosolic_introns) %>%
  left_join(N_nuclear_introns) %>%
  left_join(N_irratio) %>%
  mutate(
    cytosolic = (cyt / tot_cyt) / (high / total),
    nuclear = (nuc / tot_nuc) / (high / total)
  ) %>%
  group_by(cell_line) %>%
  summarise(
    nuc_pvalue =
      phyper((nuc - 1), high,
        low, tot_nuc,
        lower.tail = FALSE,
        log.p = FALSE
      ),
    cyt_pvalue =
      phyper((cyt - 1), high,
        low, tot_cyt,
        lower.tail = FALSE,
        log.p = FALSE
      ),
    nuclear = nuclear,
    cytosolic = cytosolic
  )

high_irratio_enr %<>%
  mutate(
    cell_line =
      str_replace_all(
        cell_line,
        "endothelial",
        "HUVEC"
      ) %>%
        str_replace_all(
          "keratinocyte",
          "NHEK"
        )
  ) %>%
  pivot_longer(4:5, names_to = "frac", values_to = "enr") %>%
  pivot_longer(2:3, names_to = "frac_pval", values_to = "pvalue") %>%
  filter((frac == "nuclear" & frac_pval == "nuc_pvalue") |
    (frac == "cytosolic" & frac_pval == "cyt_pvalue")) %>%
  select(-frac_pval) %>%
  mutate(pvalue = -log10(p.adjust(pvalue, method = "bonferroni"))) %>%
  add_column(sig_level = "") %>%
  mutate(sig_level = replace(sig_level, pvalue >= -log10(0.05), "*"))


# Plot enrichment
# ================
pdf("figures/encode_data/final/exon_skipping_intron_inclusion_enrichment/high_irratio_enrichment.pdf",
  width = 3.2,
  height = 3
)

high_irratio_enr %>%
  ggplot(aes(y = enr, x = frac, fill = cell_line)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(
    aes(
      y = enr + 0.05, x = frac,
      label = sig_level, group = cell_line
    ),
    position = position_dodge(width = .9),
    size = 2.8, family = font
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.text = element_text(size = size_medium),
    plot.tag = element_text(face = "bold"),
    axis.text = element_text(size = size_small),
    axis.text.x =
      element_text(angle = 45, vjust = 1, hjust = 1),
    legend.box.margin = margin(t = -1, b = -0.5, unit = "in"),
    plot.margin = margin(t = 1, unit = "in")
  ) +
  labs(
    y = "Enrichment of high IR-ratio introns",
    x = "Fraction specific introns",
    tag = "D"
  ) +
  scale_fill_manual(
    values = paletteer_d("ggsci::default_igv")[2:11],
    name = "ENCODE sample"
  )
dev.off()



# Plot enrichment for low psi
# and high irratio together
# ============================
pdf("figures/encode_data/final/exon_skipping_intron_inclusion_enrichment/high_irratio_low_psi_enrichment.pdf",
  width = 2,
  height = 4
)

low_psi_enr %>%
  add_column(enr_type = "Low PSI exons") %>%
  bind_rows(high_irratio_enr %>%
    add_column(enr_type = "High IR-ratio introns")) %>%
  ggplot(aes(y = enr, x = frac, fill = cell_line)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~enr_type, ncol = 1, scales = "free_y") +
  geom_text(
    aes(
      y = enr + 0.05, x = frac,
      label = sig_level, group = cell_line
    ),
    position = position_dodge(width = .9),
    size = 2.8, family = font
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(size = frame_size),
    text = element_text(size = size_big, family = font),
    legend.position = "none",
    plot.tag = element_text(face = "bold"),
    axis.text = element_text(size = size_small),
    axis.text.x =
      element_text(angle = 45, vjust = 1, hjust = 1),
    strip.text = element_text(size = size_big, face = "bold"),
    strip.background = element_rect(
      colour = "white",
      fill = "white"
    )
  ) +
  labs(
    y = "Enrichment",
    x = "Fraction specific exons/introns",
    tag = "D"
  ) +
  scale_fill_manual(values = paletteer_d("ggsci::default_igv")[2:11])
dev.off()
