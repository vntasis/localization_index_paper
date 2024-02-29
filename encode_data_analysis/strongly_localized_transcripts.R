#==========================================================
# Transcripts with strong cytosolic or nuclear localization
#==========================================================

strong_loc <-
  loc_idxs_final %>%
  mutate(`< 0.3` = ifelse(loc_idx < 0.3, TRUE, FALSE)) %>%
  mutate(`> 0.9` = ifelse(loc_idx > 0.9, TRUE, FALSE)) %>%
  pivot_longer(7:8, names_to = "cutoff", values_to = "pass") %>%
  group_by(transcript_id, cutoff) %>%
  summarise(
    transcript_type = transcript_type,
    n_cell = sum(pass)
  ) %>%
  left_join(loc_idxs_final %>%
    group_by(transcript_id) %>%
    summarise(n_expressed = n())) %>%
  distinct()

strong_cytosolic <-
  strong_loc %>%
  filter(cutoff == "> 0.9" &
    n_expressed > 4 & n_cell / n_expressed == 1) %>%
  separate(transcript_id, c("transcript_id", "x")) %>%
  ungroup() %>%
  select(transcript_id, transcript_type, n_expressed)

strong_nuclear <-
  strong_loc %>%
  filter(cutoff == "< 0.3" &
    n_expressed > 4 & n_cell / n_expressed == 1) %>%
  separate(transcript_id, c("transcript_id", "x")) %>%
  ungroup() %>%
  select(transcript_id, transcript_type, n_expressed)

write.table(strong_cytosolic,
  "highly_localized_transcripts/strong_cytosolic",
  quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = FALSE
)

write.table(strong_nuclear,
  "highly_localized_transcripts/strong_nuclear",
  quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = FALSE
)

write.table(
  strong_loc %>%
    filter(n_expressed > 4) %>%
    separate(transcript_id, c("transcript_id", "x")) %>%
    pull(transcript_id) %>%
    unique(),
  "highly_localized_transcripts/background",
  quote = FALSE, sep = "\n",
  row.names = FALSE, col.names = FALSE
)
