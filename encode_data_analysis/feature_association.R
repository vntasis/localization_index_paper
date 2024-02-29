#*********************************************************
## Descriptive stats about Localization index
## and its association with different transcript features
#*********************************************************

#---------------------
# loc_idx distribution
#---------------------

source("scripts/encode_data/locidx_distribution_per_transcript_type.R",
  echo = TRUE
)

#------------------
## Some other stats
#------------------

source("scripts/encode_data/numbers_of_nuc_cyt_transcripts.R",
  echo = TRUE
)


source("scripts/encode_data/numbers_of_consistent_nuc_cyt_transcripts_barplot.R",
  echo = TRUE
)


source("scripts/encode_data/genes_with_transcripts_of_varying_localization.R",
  echo = TRUE
)


source("scripts/encode_data/strongly_localized_transcripts.R",
  echo = TRUE
)


#---------
## Markers
#---------

source("scripts/encode_data/marker_genes_locidx.R",
  echo = TRUE
)


#------------------------------------------
## Association with transcript str features
#------------------------------------------

source("scripts/encode_data/transcript_structure_association.R",
  echo = TRUE
)


#------------------------------------------------
## Association with intron retention and splicing
## PSI and IRratio
#------------------------------------------------
source("scripts/encode_data/psi_irratio_association.R",
  echo = TRUE
)
