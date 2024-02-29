#============================================================
# How many transcripts are cytosolic in 4 cell lines or more?
#============================================================

loc_n_cell <-
  loc_idxs_final %>%
  mutate(`<= 0.5` = ifelse(loc_idx <= 0.5, TRUE, FALSE)) %>%
  mutate(`>= 0.5` = ifelse(loc_idx >= 0.5, TRUE, FALSE)) %>%
  mutate(`< 0.4` = ifelse(loc_idx < 0.4, TRUE, FALSE)) %>%
  mutate(`< 0.3` = ifelse(loc_idx < 0.3, TRUE, FALSE)) %>%
  mutate(`< 0.2` = ifelse(loc_idx < 0.2, TRUE, FALSE)) %>%
  mutate(`> 0.6` = ifelse(loc_idx > 0.6, TRUE, FALSE)) %>%
  mutate(`> 0.7` = ifelse(loc_idx > 0.7, TRUE, FALSE)) %>%
  mutate(`> 0.8` = ifelse(loc_idx > 0.8, TRUE, FALSE)) %>%
  pivot_longer(7:14, names_to = "cutoff", values_to = "pass") %>%
  group_by(transcript_id, cutoff) %>%
  summarise(n_cell = sum(pass)) %>%
  left_join(
    loc_idxs_final %>%
      group_by(transcript_id) %>%
      summarise(n_expressed = n())
  )
