# 01_QC_amino_acids.R
# Quality control of pre-filtered REGENIE mGWAS summary statistics.
# Input: one CSV per amino-acid trait in data/Amino_acids/
# Output: QC datasets and summaries in results/QC/Amino_acids/
# Note: summary-statistic QC only; no individual-level genotype QC. Amino acid as an example

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(here)

input_dir <- here("data", "Amino_acids")
output_dir <- here("results", "QC", "Amino_acids")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(input_files) == 0) stop("No input CSV files in: ", input_dir)

required_columns <- c(
  "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
  "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P",
  "EXTRA", "Pval"
)
essential_columns <- c(
  "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
  "INFO", "N", "BETA", "SE", "Pval"
)
info_threshold <- 0.80

qc_file <- function(file) {
  dat <- read_csv(file, show_col_types = FALSE)
  n_original <- nrow(dat)
  missing_columns <- setdiff(required_columns, names(dat))
  if (length(missing_columns) > 0) {
    stop(basename(file), " is missing: ", paste(missing_columns, collapse = ", "))
  }

  trait_name <- basename(file) |>
    str_remove("\\.csv$") |>
    str_remove("_variants$")

  dat <- dat |>
    mutate(
      CHROM = as.character(CHROM),
      GENPOS = as.numeric(GENPOS),
      ID = as.character(ID),
      ALLELE0 = toupper(as.character(ALLELE0)),
      ALLELE1 = toupper(as.character(ALLELE1)),
      A1FREQ = as.numeric(A1FREQ),
      INFO = as.numeric(INFO),
      N = as.numeric(N),
      TEST = as.character(TEST),
      BETA = as.numeric(BETA),
      SE = as.numeric(SE),
      CHISQ = as.numeric(CHISQ),
      LOG10P = as.numeric(LOG10P),
      Pval = as.numeric(Pval),
      trait = trait_name,
      source_file = basename(file),
      .before = 1
    ) |>
filter(if_all(all_of(essential_columns), ~ !is.na(.x))) |>
    filter(
      CHROM %in% c(as.character(1:22), "X", "Y"),
      is.finite(GENPOS), GENPOS > 0,
      str_detect(ID, "^rs[0-9]+$"),
      ALLELE0 %in% c("A", "C", "G", "T"),
      ALLELE1 %in% c("A", "C", "G", "T"),
      ALLELE0 != ALLELE1,
      is.finite(A1FREQ), between(A1FREQ, 0, 1),
      is.finite(INFO), INFO >= info_threshold,
      is.finite(N), N > 0,
      is.finite(BETA),
      is.finite(SE), SE > 0,
      is.finite(Pval), between(Pval, 0, 1),
      is.finite(LOG10P), LOG10P >= 0,
      TEST == "ADD"
    )

  # Do not silently discard conflicting duplicate trait-SNP records.
  duplicates <- dat |>
    count(trait, ID, name = "n") |>
    filter(n > 1)
  if (nrow(duplicates) > 0) {
    stop("Duplicate trait-SNP IDs found in ", basename(file),
         "; inspect them before deduplication.")
  }

  summary <- tibble(
    trait = trait_name,
    source_file = basename(file),
    n_original = n_original,
    n_after_qc = nrow(dat),
    n_removed = n_original - nrow(dat),
    proportion_retained = if (n_original > 0) nrow(dat) / n_original else NA_real_,
    n_duplicate_trait_variants = nrow(duplicates),
    info_threshold = info_threshold
  )
  list(data = dat, summary = summary)
}

qc_results <- map(input_files, qc_file)
qc_summary <- map_dfr(qc_results, "summary")
qc_data <- map(qc_results, "data")

# Save one QC file per trait, including traits with no retained records.
walk2(qc_data, qc_summary$trait, function(dat, trait) {
  write_csv(dat, file.path(output_dir, paste0(trait, "_QC.csv")))
})

amino_acids_qc <- bind_rows(qc_data) |>
  select(
    trait, CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ,
    INFO, N, TEST, BETA, SE, CHISQ, LOG10P, Pval,
    EXTRA, source_file
  )
trait_summary <- amino_acids_qc |>
  group_by(trait) |>
  summarise(
    n_variants = n(),
    n_unique_variants = n_distinct(ID),
    min_P = min(Pval),
    max_log10P = max(LOG10P),
    mean_INFO = mean(INFO),
    median_INFO = median(INFO),
    .groups = "drop"
  )

write_csv(amino_acids_qc, file.path(output_dir, "amino_acids_mGWAS_QC.csv"))
write_csv(qc_summary, file.path(output_dir, "amino_acids_QC_summary.csv"))
write_csv(trait_summary, file.path(output_dir, "amino_acids_trait_summary.csv"))

cat("QC completed.\n")
cat("Traits processed:", nrow(qc_summary), "\n")
cat("QC-passed association records:", nrow(amino_acids_qc), "\n")
cat("Unique SNP IDs:", n_distinct(amino_acids_qc$ID), "\n")
cat("Output directory:", output_dir, "\n")
print(sessionInfo())
