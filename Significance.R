# ============================================================
# 02_significant_mGWAS_associations.R
#
# Purpose:
#   Identify genome-wide significant SNP-metabolite
#   associations from QC-passed REGENIE mGWAS summary
#   statistics across six selected biochemical groups.
#
# Genome-wide significance threshold:
#   P < 5 × 10^-8
#
# Inputs:
#   QC-passed group-level mGWAS datasets
#
# Outputs:
#   1. Significant associations for each biochemical group
#   2. Combined significant association dataset
#   3. Trait-level significance summary
#   4. Biochemical-group-level significance summary
#   5. Ranked table of significant associations
#
# Important:
#   SNPs are NOT LD-clumped or reduced to independent loci
#   at this stage. Independent signal/locus evaluation will
#   be performed during the subsequent prioritisation stage.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)


# ------------------------------------------------------------
# 2. Define directories
# ------------------------------------------------------------

base_dir <- "/Users/janiyajanat/Downloads/OneDrive_1_18-08-2026"

output_dir <- file.path(
  base_dir,
  "Significant_associations"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 3. Define QC-passed input files
# ------------------------------------------------------------

input_files <- tibble(
  biochemical_group = c(
    "Amino acids",
    "Fatty acids",
    "Glycolysis",
    "Ketone bodies",
    "Lipoprotein concentrations",
    "Triglycerides"
  ),
  
  file = c(
    file.path(
      base_dir,
      "Amino_acids/QC/amino_acids_mGWAS_QC.csv"
    ),
    
    file.path(
      base_dir,
      "Fatty_acids/QC/fatty_acids_mGWAS_QC.csv"
    ),
    
    file.path(
      base_dir,
      "Glycolysis/QC/glycolysis_mGWAS_QC.csv"
    ),
    
    file.path(
      base_dir,
      "Ketone_bodies/QC/ketone_bodies_mGWAS_QC.csv"
    ),
    
    file.path(
      base_dir,
      "Lipoprotein_concentrations/QC/lipoprotein_concentrations_mGWAS_QC.csv"
    ),
    
    file.path(
      base_dir,
      "Triglycerides/QC/triglycerides_mGWAS_QC.csv"
    )
  )
)


# ------------------------------------------------------------
# 4. Confirm all files exist
# ------------------------------------------------------------

missing_files <- input_files |>
  filter(
    !file.exists(file)
  )

if (nrow(missing_files) > 0) {
  
  stop(
    paste0(
      "The following input files were not found:\n",
      paste(
        missing_files$file,
        collapse = "\n"
      )
    )
  )
}

cat(
  "All six QC-passed datasets found successfully.\n"
)


# ------------------------------------------------------------
# 5. Define genome-wide significance threshold
# ------------------------------------------------------------

GWAS_THRESHOLD <- 5e-8

cat(
  "Genome-wide significance threshold:",
  GWAS_THRESHOLD,
  "\n"
)


# ------------------------------------------------------------
# 6. Function to read each biochemical group
# ------------------------------------------------------------

read_group_data <- function(
    group_name,
    file_path
) {
  
  dat <- read_csv(
    file_path,
    show_col_types = FALSE
  )
  
  dat |>
    mutate(
      biochemical_group = group_name,
      .before = 1
    )
}


# ------------------------------------------------------------
# 7. Read all six QC-passed datasets
# ------------------------------------------------------------

all_qc_data <- map2_dfr(
  input_files$biochemical_group,
  input_files$file,
  read_group_data
)


# ------------------------------------------------------------
# 8. Confirm total input size
#
# Expected from QC:
#   44 metabolic traits
#   408,150 SNP-metabolite association records
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "INPUT DATA SUMMARY\n"
)

cat(
  "============================================\n"
)

cat(
  "Biochemical groups:",
  n_distinct(
    all_qc_data$biochemical_group
  ),
  "\n"
)

cat(
  "Metabolic traits:",
  n_distinct(
    all_qc_data$trait
  ),
  "\n"
)

cat(
  "QC-passed association records:",
  nrow(
    all_qc_data
  ),
  "\n"
)


# ------------------------------------------------------------
# 9. Check P-values before significance filtering
# ------------------------------------------------------------

pvalue_check <- all_qc_data |>
  summarise(
    missing_P = sum(
      is.na(Pval)
    ),
    
    invalid_P = sum(
      !is.na(Pval) &
        (
          Pval < 0 |
            Pval > 1
        )
    ),
    
    minimum_P = min(
      Pval,
      na.rm = TRUE
    ),
    
    maximum_P = max(
      Pval,
      na.rm = TRUE
    )
  )

cat(
  "\nP-value check:\n"
)

print(
  pvalue_check
)


# ------------------------------------------------------------
# 10. Identify genome-wide significant associations
#
# IMPORTANT:
#   This identifies significant SNP-metabolite associations.
#   It does NOT yet identify independent loci.
# ------------------------------------------------------------

significant_associations <- all_qc_data |>
  filter(
    Pval < GWAS_THRESHOLD
  )


# ------------------------------------------------------------
# 11. Rank all significant associations
#
# Smallest P-value = strongest statistical association
# ------------------------------------------------------------

significant_associations <- significant_associations |>
  arrange(
    Pval
  )


# ------------------------------------------------------------
# 12. Save complete significant-association dataset
# ------------------------------------------------------------

write_csv(
  significant_associations,
  file.path(
    output_dir,
    "all_groups_mGWAS_significant.csv"
  )
)


# ------------------------------------------------------------
# 13. Save significant associations separately by group
# ------------------------------------------------------------

group_file_names <- c(
  "Amino acids" =
    "amino_acids_mGWAS_significant.csv",
  
  "Fatty acids" =
    "fatty_acids_mGWAS_significant.csv",
  
  "Glycolysis" =
    "glycolysis_mGWAS_significant.csv",
  
  "Ketone bodies" =
    "ketone_bodies_mGWAS_significant.csv",
  
  "Lipoprotein concentrations" =
    "lipoprotein_concentrations_mGWAS_significant.csv",
  
  "Triglycerides" =
    "triglycerides_mGWAS_significant.csv"
)


for (
  group_name in names(group_file_names)
) {
  
  group_data <- significant_associations |>
    filter(
      biochemical_group == group_name
    )
  
  write_csv(
    group_data,
    file.path(
      output_dir,
      group_file_names[[group_name]]
    )
  )
}


# ------------------------------------------------------------
# 14. Generate trait-level summary
#
# IMPORTANT:
#   n_significant_snps represents significant association
#   records, NOT independent genetic loci.
# ------------------------------------------------------------

trait_summary <- all_qc_data |>
  group_by(
    biochemical_group,
    trait
  ) |>
  summarise(
    
    n_tested_after_qc = n(),
    
    n_significant_snps = sum(
      Pval < GWAS_THRESHOLD,
      na.rm = TRUE
    ),
    
    strongest_P = min(
      Pval,
      na.rm = TRUE
    ),
    
    max_log10P = max(
      LOG10P,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  
  mutate(
    
    has_genome_wide_signal =
      n_significant_snps > 0,
    
    percentage_significant =
      100 *
      n_significant_snps /
      n_tested_after_qc
  ) |>
  
  arrange(
    desc(n_significant_snps),
    strongest_P
  )


# ------------------------------------------------------------
# 15. Save trait-level summary
# ------------------------------------------------------------

write_csv(
  trait_summary,
  file.path(
    output_dir,
    "trait_level_significance_summary.csv"
  )
)


# ------------------------------------------------------------
# 16. Generate group-level summary
# ------------------------------------------------------------

group_summary <- all_qc_data |>
  group_by(
    biochemical_group
  ) |>
  summarise(
    
    n_traits =
      n_distinct(trait),
    
    n_qc_passed_associations =
      n(),
    
    n_significant_associations =
      sum(
        Pval < GWAS_THRESHOLD,
        na.rm = TRUE
      ),
    
    n_traits_with_significant_signal =
      n_distinct(
        trait[
          Pval < GWAS_THRESHOLD
        ]
      ),
    
    strongest_P =
      min(
        Pval,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  
  mutate(
    
    percentage_significant =
      100 *
      n_significant_associations /
      n_qc_passed_associations
  ) |>
  
  arrange(
    desc(
      n_significant_associations
    )
  )


# ------------------------------------------------------------
# 17. Save group-level summary
# ------------------------------------------------------------

write_csv(
  group_summary,
  file.path(
    output_dir,
    "biochemical_group_significance_summary.csv"
  )
)


# ------------------------------------------------------------
# 18. Generate ranked table of significant associations
# ------------------------------------------------------------

ranked_significant <- significant_associations |>
  select(
    biochemical_group,
    trait,
    CHROM,
    GENPOS,
    ID,
    ALLELE0,
    ALLELE1,
    A1FREQ,
    INFO,
    N,
    BETA,
    SE,
    Pval,
    LOG10P
  ) |>
  arrange(
    Pval
  ) |>
  mutate(
    overall_rank = row_number(),
    .before = 1
  )


# ------------------------------------------------------------
# 19. Save ranked significant associations
# ------------------------------------------------------------

write_csv(
  ranked_significant,
  file.path(
    output_dir,
    "ranked_significant_associations.csv"
  )
)


# ------------------------------------------------------------
# 20. Generate top association for each metabolic trait
# ------------------------------------------------------------

top_hit_per_trait <- significant_associations |>
  group_by(
    biochemical_group,
    trait
  ) |>
  slice_min(
    order_by = Pval,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  arrange(
    Pval
  )


# ------------------------------------------------------------
# 21. Save top hit per trait
# ------------------------------------------------------------

write_csv(
  top_hit_per_trait,
  file.path(
    output_dir,
    "top_significant_hit_per_trait.csv"
  )
)


# ------------------------------------------------------------
# 22. Generate summary of significant traits only
# ------------------------------------------------------------

significant_trait_summary <- trait_summary |>
  filter(
    n_significant_snps > 0
  ) |>
  arrange(
    desc(n_significant_snps),
    strongest_P
  )

write_csv(
  significant_trait_summary,
  file.path(
    output_dir,
    "significant_traits_only_summary.csv"
  )
)


# ------------------------------------------------------------
# 23. Final significance summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "GENOME-WIDE SIGNIFICANCE SUMMARY\n"
)

cat(
  "============================================\n"
)

cat(
  "Genome-wide threshold: P <",
  GWAS_THRESHOLD,
  "\n"
)

cat(
  "Total QC-passed associations tested:",
  nrow(
    all_qc_data
  ),
  "\n"
)

cat(
  "Total genome-wide significant associations:",
  nrow(
    significant_associations
  ),
  "\n"
)

cat(
  "Metabolic traits analysed:",
  n_distinct(
    all_qc_data$trait
  ),
  "\n"
)

cat(
  "Traits with >=1 genome-wide significant association:",
  n_distinct(
    significant_associations$trait
  ),
  "\n"
)

cat(
  "Traits without genome-wide significant associations:",
  n_distinct(
    all_qc_data$trait
  ) -
    n_distinct(
      significant_associations$trait
    ),
  "\n"
)

cat(
  "Unique SNPs among significant associations:",
  n_distinct(
    significant_associations$ID
  ),
  "\n"
)


# ------------------------------------------------------------
# 24. Display biochemical-group summary
# ------------------------------------------------------------

cat(
  "\nSignificant associations by biochemical group:\n"
)

print(
  group_summary,
  n = Inf
)


# ------------------------------------------------------------
# 25. Display trait-level summary
# ------------------------------------------------------------

cat(
  "\nSignificance summary by metabolic trait:\n"
)

print(
  trait_summary,
  n = Inf
)


# ------------------------------------------------------------
# 26. Display strongest associations
# ------------------------------------------------------------

cat(
  "\nTop 20 genome-wide significant associations:\n"
)

print(
  ranked_significant |>
    slice_head(
      n = 20
    )
)


# ------------------------------------------------------------
# 27. End of script
# ------------------------------------------------------------

cat(
  "\nSignificance filtering completed successfully.\n"
)

cat(
  "Results saved to:\n",
  output_dir,
  "\n"
)