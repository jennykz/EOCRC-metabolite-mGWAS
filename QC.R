# ============================================================
# 01_QC_amino_acids.R
#
# Purpose:
#   Quality control and standardisation of REGENIE mGWAS
#   summary statistics for amino-acid metabolic traits.
#
# Input:
#   One CSV file per metabolic trait.
#
# Output:
#   1. QC'd individual trait files
#   2. Combined QC'd amino-acid dataset
#   3. QC summary table
#   4. Trait-level variant summary table
#
# QC performed:
#   - File and column structure
#   - Missing essential values
#   - Chromosome and genomic position
#   - Variant ID format
#   - SNP allele format
#   - Allele frequency range
#   - INFO score
#   - Association statistics
#   - Additive association model
#   - Duplicate trait-SNP combinations
#
# Note:
#   This script performs QC of mGWAS summary statistics
#   generated upstream. It does NOT perform QC of individual-
#   level UK Biobank data.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)


# ------------------------------------------------------------
# 2. Define input and output directories
# ------------------------------------------------------------

input_dir <- "/Users/janiyajanat/Downloads/OneDrive_1_18-08-2026/Amino_acids"

output_dir <- file.path(
  input_dir,
  "QC"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Identify raw mGWAS CSV files
# ------------------------------------------------------------

input_files <- list.files(
  path = input_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Exclude CSV files located in the QC directory
input_files <- input_files[
  !str_detect(input_files, "/QC/")
]

cat(
  "Number of input files:",
  length(input_files),
  "\n"
)

print(
  basename(input_files)
)


# ------------------------------------------------------------
# 4. Define expected columns
# ------------------------------------------------------------

expected_columns <- c(
  "CHROM",
  "GENPOS",
  "ID",
  "ALLELE0",
  "ALLELE1",
  "A1FREQ",
  "INFO",
  "N",
  "TEST",
  "BETA",
  "SE",
  "CHISQ",
  "LOG10P",
  "EXTRA",
  "Pval"
)


# ------------------------------------------------------------
# 5. Define essential columns
#
# These variables must contain valid values for an
# association record to be retained.
# ------------------------------------------------------------

essential_columns <- c(
  "CHROM",
  "GENPOS",
  "ID",
  "ALLELE0",
  "ALLELE1",
  "INFO",
  "N",
  "BETA",
  "SE",
  "Pval"
)


# ------------------------------------------------------------
# 6. Define QC thresholds
# ------------------------------------------------------------

INFO_THRESHOLD <- 0.80


# ------------------------------------------------------------
# 7. Define QC function
# ------------------------------------------------------------

qc_mgwas_file <- function(
    file,
    info_threshold = INFO_THRESHOLD
) {
  
  # ----------------------------------------------------------
  # Read mGWAS summary statistics
  # ----------------------------------------------------------
  
  dat <- read_csv(
    file,
    show_col_types = FALSE
  )
  
  n_original <- nrow(dat)
  
  
  # ----------------------------------------------------------
  # Check expected columns
  # ----------------------------------------------------------
  
  missing_columns <- setdiff(
    expected_columns,
    colnames(dat)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste0(
        "\nFile: ",
        basename(file),
        "\nMissing columns: ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        "\n"
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Extract metabolic trait name from filename
  #
  # Example:
  # assoc.23460_alanine_variants.csv
  #
  # becomes:
  # assoc.23460_alanine
  # ----------------------------------------------------------
  
  trait_name <- basename(file) |>
    str_remove("\\.csv$") |>
    str_remove("_variants$")
  
  
  # ----------------------------------------------------------
  # Standardise variable types
  # ----------------------------------------------------------
  
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
      Pval = as.numeric(Pval)
    )
  
  
  # ----------------------------------------------------------
  # Add metabolic trait and source-file identifiers
  # ----------------------------------------------------------
  
  dat <- dat |>
    mutate(
      trait = trait_name,
      source_file = basename(file),
      .before = 1
    )
  
  
  # ----------------------------------------------------------
  # QC 1: Remove records with missing essential fields
  # ----------------------------------------------------------
  
  dat <- dat |>
    filter(
      if_all(
        all_of(essential_columns),
        ~ !is.na(.x)
      )
    )
  
  
  # ----------------------------------------------------------
  # QC 2: Check chromosome and genomic position
  # ----------------------------------------------------------
  
  valid_chromosomes <- c(
    as.character(1:22),
    "X",
    "Y"
  )
  
  dat <- dat |>
    filter(
      CHROM %in% valid_chromosomes,
      GENPOS > 0
    )
  
  
  # ----------------------------------------------------------
  # QC 3: Retain valid dbSNP rsIDs
  # ----------------------------------------------------------
  
  dat <- dat |>
    filter(
      str_detect(
        ID,
        "^rs[0-9]+$"
      )
    )
  
  
  # ----------------------------------------------------------
  # QC 4: Retain biallelic SNPs
  #
  # Only single-nucleotide A/C/G/T alleles are retained.
  # ----------------------------------------------------------
  
  valid_bases <- c(
    "A",
    "C",
    "G",
    "T"
  )
  
  dat <- dat |>
    filter(
      ALLELE0 %in% valid_bases,
      ALLELE1 %in% valid_bases,
      nchar(ALLELE0) == 1,
      nchar(ALLELE1) == 1,
      ALLELE0 != ALLELE1
    )
  
  
  # ----------------------------------------------------------
  # QC 5: Validate allele frequency
  #
  # No additional MAF threshold is applied.
  # Values are checked only to ensure they fall within
  # the valid frequency range of 0-1.
  # ----------------------------------------------------------
  
  dat <- dat |>
    filter(
      is.finite(A1FREQ),
      A1FREQ >= 0,
      A1FREQ <= 1
    )
  
  
  # ----------------------------------------------------------
  # QC 6: Filter variants according to INFO score
  #
  # Variants with INFO < 0.80 are excluded.
  # ----------------------------------------------------------
  
  dat <- dat |>
    filter(
      is.finite(INFO),
      INFO >= info_threshold
    )
  
  
  # ----------------------------------------------------------
  # QC 7: Validate association statistics
  # ----------------------------------------------------------
  
  dat <- dat |>
    filter(
      is.finite(BETA),
      is.finite(SE),
      SE > 0,
      is.finite(Pval),
      Pval >= 0,
      Pval <= 1,
      is.finite(LOG10P),
      LOG10P >= 0
    )
  
  
  # ----------------------------------------------------------
  # QC 8: Retain additive association tests
  #
  # The upstream REGENIE results use the ADD test.
  # ----------------------------------------------------------
  
  dat <- dat |>
    filter(
      TEST == "ADD"
    )
  
  
  # ----------------------------------------------------------
  # QC 9: Identify duplicate trait-SNP combinations
  #
  # The same SNP may legitimately occur for multiple
  # metabolic traits. Therefore, SNP ID alone is not used
  # to remove duplicates.
  # ----------------------------------------------------------
  
  duplicate_trait_variants <- dat |>
    count(
      trait,
      ID,
      name = "n"
    ) |>
    filter(
      n > 1
    )
  
  
  # Retain one record for each trait-SNP combination
  dat <- dat |>
    distinct(
      trait,
      ID,
      .keep_all = TRUE
    )
  
  
  # ----------------------------------------------------------
  # Generate QC summary for this metabolic trait
  # ----------------------------------------------------------
  
  qc_summary <- tibble(
    trait = trait_name,
    source_file = basename(file),
    n_original = n_original,
    n_after_qc = nrow(dat),
    n_removed = n_original - nrow(dat),
    proportion_retained =
      nrow(dat) / n_original,
    n_duplicate_trait_variants =
      nrow(duplicate_trait_variants),
    info_threshold = info_threshold
  )
  
  
  # ----------------------------------------------------------
  # Return QC'd data and summary
  # ----------------------------------------------------------
  
  list(
    data = dat,
    summary = qc_summary
  )
}


# ------------------------------------------------------------
# 8. Run QC on all amino-acid files
# ------------------------------------------------------------

qc_results <- map(
  input_files,
  qc_mgwas_file
)


# ------------------------------------------------------------
# 9. Extract QC'd datasets
# ------------------------------------------------------------

qc_data <- map(
  qc_results,
  "data"
)


# ------------------------------------------------------------
# 10. Extract QC summaries
# ------------------------------------------------------------

qc_summary <- map_dfr(
  qc_results,
  "summary"
)


# ------------------------------------------------------------
# 11. Save individual QC'd trait files
# ------------------------------------------------------------

walk(
  qc_data,
  function(dat) {
    
    trait_name <- unique(
      dat$trait
    )
    
    output_file <- file.path(
      output_dir,
      paste0(
        trait_name,
        "_QC.csv"
      )
    )
    
    write_csv(
      dat,
      output_file
    )
  }
)


# ------------------------------------------------------------
# 12. Combine all QC'd amino-acid traits
# ------------------------------------------------------------

amino_acids_QC <- bind_rows(
  qc_data
)


# ------------------------------------------------------------
# 13. Organise columns
# ------------------------------------------------------------

amino_acids_QC <- amino_acids_QC |>
  select(
    trait,
    CHROM,
    GENPOS,
    ID,
    ALLELE0,
    ALLELE1,
    A1FREQ,
    INFO,
    N,
    TEST,
    BETA,
    SE,
    CHISQ,
    LOG10P,
    Pval,
    EXTRA,
    source_file
  )


# ------------------------------------------------------------
# 14. Save combined QC'd dataset
# ------------------------------------------------------------

write_csv(
  amino_acids_QC,
  file.path(
    output_dir,
    "amino_acids_mGWAS_QC.csv"
  )
)


# ------------------------------------------------------------
# 15. Save QC summary
# ------------------------------------------------------------

write_csv(
  qc_summary,
  file.path(
    output_dir,
    "amino_acids_QC_summary.csv"
  )
)


# ------------------------------------------------------------
# 16. Generate trait-level variant summary
# ------------------------------------------------------------

trait_variant_counts <- amino_acids_QC |>
  group_by(
    trait
  ) |>
  summarise(
    n_variants = n(),
    n_unique_variants = n_distinct(ID),
    min_P = min(
      Pval,
      na.rm = TRUE
    ),
    max_log10P = max(
      LOG10P,
      na.rm = TRUE
    ),
    mean_INFO = mean(
      INFO,
      na.rm = TRUE
    ),
    median_INFO = median(
      INFO,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(
    desc(n_variants)
  )


# ------------------------------------------------------------
# 17. Save trait-level summary
# ------------------------------------------------------------

write_csv(
  trait_variant_counts,
  file.path(
    output_dir,
    "amino_acids_trait_summary.csv"
  )
)


# ------------------------------------------------------------
# 18. Final QC checks
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "FINAL QC SUMMARY\n"
)

cat(
  "============================================\n"
)


# Number of metabolic traits
cat(
  "Number of metabolic traits:",
  n_distinct(
    amino_acids_QC$trait
  ),
  "\n"
)


# Number of SNP-metabolite association records
cat(
  "Total QC-passed SNP-metabolite associations:",
  nrow(amino_acids_QC),
  "\n"
)


# Number of unique SNPs
cat(
  "Unique SNPs across all traits:",
  n_distinct(
    amino_acids_QC$ID
  ),
  "\n"
)


# ------------------------------------------------------------
# 19. Check for missing values in essential fields
# ------------------------------------------------------------

missing_check <- amino_acids_QC |>
  summarise(
    across(
      all_of(essential_columns),
      ~ sum(is.na(.x))
    )
  )

cat(
  "\nMissing values in essential fields:\n"
)

print(
  missing_check
)


# ------------------------------------------------------------
# 20. Check for duplicate trait-SNP combinations
# ------------------------------------------------------------

duplicate_check <- amino_acids_QC |>
  count(
    trait,
    ID
  ) |>
  filter(
    n > 1
  )

cat(
  "\nDuplicate trait-SNP combinations:",
  nrow(duplicate_check),
  "\n"
)


# ------------------------------------------------------------
# 21. Check sample size consistency
# ------------------------------------------------------------

sample_size_check <- amino_acids_QC |>
  count(
    trait,
    N
  ) |>
  arrange(
    trait,
    N
  )

cat(
  "\nSample size by trait:\n"
)

print(
  sample_size_check
)


# ------------------------------------------------------------
# 22. Summarise INFO scores
# ------------------------------------------------------------

info_summary <- amino_acids_QC |>
  summarise(
    minimum_INFO = min(
      INFO,
      na.rm = TRUE
    ),
    median_INFO = median(
      INFO,
      na.rm = TRUE
    ),
    mean_INFO = mean(
      INFO,
      na.rm = TRUE
    ),
    maximum_INFO = max(
      INFO,
      na.rm = TRUE
    )
  )

cat(
  "\nINFO score summary:\n"
)

print(
  info_summary
)


# ------------------------------------------------------------
# 23. Display trait-level summary
# ------------------------------------------------------------

cat(
  "\nTrait-level variant summary:\n"
)

print(
  trait_variant_counts
)


# ------------------------------------------------------------
# End of script
# ------------------------------------------------------------
