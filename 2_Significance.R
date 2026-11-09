# 02_significant_mGWAS_associations.R
# Identify genome-wide significant SNP-metabolite associations (P < 5e-8).
# Input: QC-passed group datasets in results/QC/<group>/
# Output: significant associations and summaries in results/significant_associations/
# Upstream REGENIE summary statistics were pre-filtered at P <= 1e-3.
# Significant association records are NOT independent LD loci.

library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(here)

threshold <- 5e-8
output_dir <- here("results", "significant_associations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

groups <- tibble(
  biochemical_group = c(
    "Amino acids", "Fatty acids", "Glycolysis", "Ketone bodies",
    "Lipoprotein concentrations", "Triglycerides"
  ),
  folder = c(
    "Amino_acids", "Fatty_acids", "Glycolysis", "Ketone_bodies",
    "Lipoprotein_concentrations", "Triglycerides"
  ),
  prefix = c(
    "amino_acids", "fatty_acids", "glycolysis", "ketone_bodies",
    "lipoprotein_concentrations", "triglycerides"
  )
) |>
  mutate(file = map2_chr(folder, prefix, ~ here(
    "results", "QC", .x, paste0(.y, "_mGWAS_QC.csv")
  )))

missing_files <- groups$file[!file.exists(groups$file)]
if (length(missing_files) > 0) {
  stop("Missing QC files:\n", paste(missing_files, collapse = "\n"))
}

read_group <- function(file, group) {
  dat <- read_csv(file, show_col_types = FALSE)
  required <- c("trait", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
                "A1FREQ", "INFO", "N", "BETA", "SE", "Pval", "LOG10P")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop(basename(file), " is missing: ", paste(missing, collapse = ", "))
  }
  if (any(!is.finite(dat$Pval) | dat$Pval < 0 | dat$Pval > 1)) {
    stop("Invalid or missing P-values in ", basename(file))
  }
  mutate(dat, biochemical_group = group, .before = 1)
}

all_qc <- map2_dfr(groups$file, groups$biochemical_group, read_group)
if (nrow(all_qc) == 0) stop("No QC-passed associations available.")

significant <- all_qc |>
  filter(Pval < threshold) |>
  arrange(Pval)

write_csv(significant, file.path(output_dir, "all_groups_mGWAS_significant.csv"))

# Save one significant-association dataset per biochemical group.
pwalk(groups[c("biochemical_group", "prefix")], function(biochemical_group, prefix) {
  significant |>
    filter(.data$biochemical_group == .env$biochemical_group) |>
    write_csv(file.path(output_dir, paste0(prefix, "_mGWAS_significant.csv")))
})

trait_summary <- all_qc |>
  group_by(biochemical_group, trait) |>
  summarise(
    n_available_after_qc = n(),
    n_significant_snps = sum(Pval < threshold),
    strongest_P = min(Pval),
    max_log10P = max(LOG10P),
    .groups = "drop"
  ) |>
  mutate(
    has_genome_wide_signal = n_significant_snps > 0,
    percentage_significant_in_available_data =
      100 * n_significant_snps / n_available_after_qc
  ) |>
  arrange(desc(n_significant_snps), strongest_P)

group_summary <- all_qc |>
  group_by(biochemical_group) |>
  summarise(
    n_traits = n_distinct(trait),
    n_qc_passed_associations = n(),
    n_significant_associations = sum(Pval < threshold),
    n_traits_with_significant_signal = n_distinct(trait[Pval < threshold]),
    strongest_P = min(Pval),
    .groups = "drop"
  ) |>
  mutate(percentage_significant_in_available_data =
           100 * n_significant_associations / n_qc_passed_associations) |>
  arrange(desc(n_significant_associations))

ranked <- significant |>
  select(biochemical_group, trait, CHROM, GENPOS, ID, ALLELE0,
         ALLELE1, A1FREQ, INFO, N, BETA, SE, Pval, LOG10P) |>
  mutate(overall_rank = row_number(), .before = 1)

top_per_trait <- significant |>
  group_by(biochemical_group, trait) |>
  slice_min(Pval, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(Pval)

write_csv(trait_summary, file.path(output_dir, "trait_level_significance_summary.csv"))
write_csv(group_summary, file.path(output_dir, "biochemical_group_significance_summary.csv"))
write_csv(ranked, file.path(output_dir, "ranked_significant_associations.csv"))
write_csv(top_per_trait, file.path(output_dir, "top_significant_hit_per_trait.csv"))
write_csv(filter(trait_summary, n_significant_snps > 0),
          file.path(output_dir, "significant_traits_only_summary.csv"))

cat("QC-passed association records:", nrow(all_qc), "\n")
cat("Genome-wide significant associations:", nrow(significant), "\n")
cat("Traits analysed:", n_distinct(all_qc$trait), "\n")
cat("Traits with significant associations:", n_distinct(significant$trait), "\n")
cat("Output directory:", output_dir, "\n")
print(sessionInfo())
