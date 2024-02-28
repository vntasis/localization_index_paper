#!/usr/bin/env r

#****************************************************
# This script takes a vector of integer numbers
# and transforms it into two vectors of equal length,
# whose sum equals the original vector.
#****************************************************

#---------------
## Initial setup
#---------------

### Load the required libraries
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(vroom)))

### Disable scientific notation
options(scipen = 999)

### Argument input
if (length(argv) == 0) {
  # default arguments
  split <- 0.5
  annotation_file <- "annotation_transcripts.tsv"
} else if (length(argv) == 2) {
  split <- as.numeric(argv[1])
  annotation_file <- argv[2]
} else {
  stop("Wrong arguments.", call. = FALSE)
}

#----------------------
## Load annotation file
#----------------------

annotation <- vroom(annotation_file)
annotation %<>%
  select(transcript_id, total_exon_len)

## ------------------------------------
## Define the transformation function.
## ------------------------------------

split2sum <-
  function(N, spl = 0.5) {
    # Divides a number into two parts.
    # Uses negative binomial to pick the difference between the two.
    # A + B = N
    # A - B = K

    out <- integer()
    a_m <- round(N * split)
    k_m <- abs(2 * a_m - N)
    probs <- dnbinom(0:N, 10, mu = k_m)
    probs <- probs / sum(probs)
    k <- sample(0:N, 1, prob = probs)
    if (split > 0.5) {
      out[1] <- round((N + k) / 2)
    } else {
      out[2] <- round((N + k) / 2)
    }
    if (split > 0.5) out[2] <- N - out[1] else out[1] <- N - out[2]
    out
  }

## ------------------
## Split the numbers
## ------------------

n_wcell <- X$N

set.seed(123)

print("Spliting...")
x <-
  sapply(n_wcell, split2sum, spl = split) %>%
  t() %>%
  as_tibble()



### Create the relative frequencies (proportions) per transcript

x %<>%
  mutate(
    pV1 = V1 / sum(V1),
    pV2 = V2 / sum(V2),
    split = V2 / (V1 + V2)
  )



### Save the generated numbers

print("Saving...")
write.table(x[, c(3, 1)], "n_nuclear.tsv", quote = FALSE,
            sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(x[, c(4, 2)], "n_cytosolic.tsv", quote = FALSE,
            sep = "\t", col.names = FALSE, row.names = FALSE)


#----------------------------------------------------
## Calculate beta and print summary of the simulation
#----------------------------------------------------

beta <-
  X %>%
  left_join(annotation) %>%
  cbind(x) %>%
  as_tibble() %>%
  mutate(num = total_exon_len * V2, denom = total_exon_len * N) %$%
  `/`(sum(num), sum(denom))

print(paste0("simulated beta: ", round(beta, 3)))

print(
  x %>%
    # drop_na() %>%
    pull(split) %>%
    summary()
)


#--------------------------------------------------------------
## Make a histogram as an illustrated summary of the simulation
#--------------------------------------------------------------

pdf("loc_idx_histogram.pdf")
opar <- par(lwd = 3.5)
hist_plot <-
  x %>%
  pull(split) %>%
  hist(main = "", col = "#a6bddb", border = "black")

label_pos <-
  hist_plot %$%
  counts %>%
  max() %>%
  `*`(0.8)

abline(v = beta, col = "red", lwd = 2)
text(x = beta + 0.1, y = label_pos,
     labels = paste0("beta: ", round(beta, 3)), col = "red")
dev.off()
