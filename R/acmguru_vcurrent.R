# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(grid)

# Function to read data files and merge with sample phenotype information
read_data_file <- function(input_path, samples_file_path) {
  files <- character()

  if (dir.exists(input_path)) {
    files <-
      list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
  } else if (file.exists(input_path) &&
             grepl("\\.csv$", input_path)) {
    files <- input_path
  } else {
    stop("Input path is neither a valid file nor a directory.")
  }

  cat("Reading samples information from: ", samples_file_path, "\n")
  samples_df <-
    read.table(
      samples_file_path,
      sep = "\t",
      col.names = c("sample", "cohort_pheno"),
      header = FALSE
    )

  data_list <- list()
  for (file_path in files) {
    cat("Reading file: ", file_path, "\n")
    df <- read.csv(file_path, sep = ",")
    df <- merge(df, samples_df, by = "sample", all.x = TRUE)
    data_list[[basename(file_path)]] <- df
  }

  return(data_list)
}

# comp_het_flag ----
# Function to set compound heterozygous flags.  WARNING NOT PHASE CHECKED
set_comp_het_flag <- function(df) {
  cat("Setting compound heterozygous flag...\n")
  df <- df %>%
    group_by(sample, SYMBOL) %>%
    mutate(comp_het_flag = ifelse(n() > 1, 1, NA)) %>%
    ungroup() %>%
    mutate(comp_het_flag = ifelse(is.na(comp_het_flag) &
                                    genotype == 2, 1, comp_het_flag))

  return(df)
}

# Function to filter data based on allele frequency
filter_by_allele_frequency <- function(df, af_threshold) {
  df_filtered <- df %>%
    dplyr::filter(AF.x < af_threshold)

  return(df_filtered)
}

# PVS1 ----
# PVS1 are null variants where IMPACT=="HIGH" and inheritance match, in gene where LoF cause disease.
# Function to apply ACMG PVS1 criterion
apply_acmg_pvs1 <- function(df) {
  cat("Applying ACMG PVS1 criterion...\n")

  df <- df %>%
    mutate(ACMG_PVS1 = NA)

  # Always apply this step as it does not depend on 'Inheritance'
  df <- df %>%
    mutate(ACMG_PVS1 = ifelse(IMPACT == "HIGH" & genotype == 2, "PVS1", ACMG_PVS1))

  # Apply 'Inheritance' related step only if the column exists
  if ("Inheritance" %in% colnames(df)) {
    df <- df %>%
      mutate(ACMG_PVS1 = ifelse(IMPACT == "HIGH" & Inheritance == "AD", "PVS1", ACMG_PVS1))
  }

  # Apply comp_het related step
  df <- df %>%
    group_by(sample, SYMBOL) %>%
    mutate(ACMG_PVS1 = ifelse(ACMG_PVS1 == "PVS1", "PVS1",
                              ifelse(sum(IMPACT == "HIGH" & comp_het_flag == 1, na.rm = TRUE) >= 2 & IMPACT == "HIGH", "PVS1", ACMG_PVS1))) %>%
    ungroup()

  return(df)
}

# PS1 ----
# PS1 Same amino acid change as a previously established pathogenic variant regardless of nucleotide change. Note to keep splicing variant as PSV1 (these are covered by IPACT HIGH).
apply_acmg_ps1 <- function(df) {
  cat("Applying ACMG PS1 criterion...\n")

  if ("CLIN_SIG" %in% colnames(df)) {
    df <- df %>%
      mutate(ACMG_PS1 = ifelse(CLIN_SIG == "pathogenic", "PS1", NA))
  } else {
    cat("Column 'CLIN_SIG' not found. Skipping PS1.\n")
  }

  return(df)
}


# PS2 skip ----
# PS2 De novo (both maternity and paternity confirmed) in a patient with the disease and no family history
# Skip due to no parental genetics.

# PS3 skip ----
# PS3 Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.
# df$ACMG <- ifelse(uniprot == "pathogenic" |
# 							pubmed == "pathogenic",
# 						"PS3")

# PS4 skip ----
# The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls
# Skip do to statistical analysis separately

# PS5 ----
# The user has additional (value) strong pathogenic evidence
apply_acmg_ps5 <- function(df) {
  cat("Applying ACMG PS5 criterion...\n")

  # Check if the necessary columns exist
  if ("IMPACT" %in% colnames(df)) {
    df <- df %>%
      group_by(sample, SYMBOL) %>%
      mutate(ACMG_PS5 = ifelse(any(IMPACT == "HIGH") & comp_het_flag == 1, "PS5", NA)) %>%
      ungroup()
  } else {
    cat("Necessary columns for PS5 not found. Skipping PS5.\n")
  }

  return(df)
}


# PM2 ----
# Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium

apply_acmg_pm2 <- function(df) {
  cat("Applying ACMG PM2 criterion...\n")

  if ("gnomAD_AF" %in% colnames(df)) {
    df$gnomAD_AF <- as.numeric(df$gnomAD_AF)
    # gnomad_max <- 1e-6  # Define the maximum frequency considered for PM2
    df <- df %>%
      mutate(ACMG_PM2 = ifelse(gnomAD_AF < gnomad_max, "PM2", NA))
  } else {
    cat("Column 'gnomAD_AF' not found. Skipping PM2.\n")
  }

  return(df)
}


# PM3 ----
# For recessive disorders, detected in trans with a pathogenic variant
# some redundancy with our PS5 since our rare disease cohort filtering call IMPACT==HIGH equates pathogenic

apply_acmg_pm3 <- function(df) {
  cat("Applying ACMG PM3 criterion...\n")

  # Assuming comp_het_flag, ACMG_PS1, and ACMG_PS5 are already calculated and present in df
  if (all(c("comp_het_flag", "ACMG_PS1", "ACMG_PS5") %in% colnames(df))) {
    df <- df %>%
      group_by(sample, SYMBOL) %>%
      mutate(ACMG_PM3 = ifelse(comp_het_flag == 1 & (ACMG_PS1 == "PS1" | ACMG_PS5 == "PS5"), "PM3", NA)) %>%
      ungroup()
  } else {
    cat("Necessary columns for PM3 not found. Skipping PM3.\n")
  }

  return(df)
}

# PP3 ----
# Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.). CUSTOM: assigned if >=3 thresholds passed.
# In-Silico Predictions: VarSome now implements the ClinGen recommendations from Evidence-based calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for clinical use of PP3/BP4 criteria: # Only one engine at a time is used, depending on availability of data, in order: MitoTip & MitImpact, MetaRNN, CADD (Premium only), DANN (if CADD is not available). # The maximum strength allowed for rules PP3 & BP4 is Strong, even if there may be evidence for Very Strong, with the exception of variants that are predicted splicing (ie: similar to PVS1). # The strength is limited to Supporting, if there's Moderate evidence from rules PM1 or PM5. # Splice prediction (scSNV) is given priority over the other in-silico predictions. # conservation is used for some low-sensitivity variant types, or if no other in-silico prediction is available. Please refer to PP3 and BP4 for more specific detail.

apply_acmg_pp3 <- function(df) {
  cat("Applying ACMG PP3 criterion...\n")

  # Check if the necessary columns exist
  required_columns <- c("CADD_PHRED", "REVEL_rankscore", "MetaLR_pred",
                        "MutationAssessor_pred", "SIFT", "PolyPhen")
  if(!all(required_columns %in% colnames(df))) {
    cat("One or more required columns for PP3 are missing. Skipping PP3.\n")
    return(df)
  }

  # Prepare data for conditions
  df <- df %>%
    separate(SIFT, into = c("SIFT_label", "SIFT_score"), sep = "\\(", remove = TRUE) %>%
    mutate(SIFT_score = str_replace(SIFT_score, "\\)", "")) %>%
    separate(PolyPhen, into = c("PolyPhen_label", "PolyPhen_score"), sep = "\\(", remove = TRUE) %>%
    mutate(PolyPhen_score = str_replace(PolyPhen_score, "\\)", "")) %>%
    mutate(
      CADD_PHRED = as.numeric(CADD_PHRED),
      REVEL_rankscore = as.numeric(REVEL_rankscore),
      ACMG_PP3 = NA
    )

  # Define conditions
  cond_CADD_PHRED <- df$CADD_PHRED >= 30
  cond_REVEL_rankscore <- df$REVEL_rankscore > 0.5
  cond_MetaLR_pred <- df$MetaLR_pred == "D"
  cond_MutationAssessor_pred <- df$MutationAssessor_pred == "H"
  cond_SIFT_label <- df$SIFT_label == "deleterious"
  cond_PolyPhen_label <- df$PolyPhen_label == "probably_damaging"

  # Count the conditions met
  df$ACMG_PP3_count <- rowSums(cbind(cond_CADD_PHRED, cond_REVEL_rankscore,
                                     cond_MetaLR_pred, cond_MutationAssessor_pred,
                                     cond_SIFT_label, cond_PolyPhen_label), na.rm = TRUE)

  # Apply PP3 if the threshold is met
  threshold <- 3
  df$ACMG_PP3 <- ifelse(df$ACMG_PP3_count >= threshold, "PP3", NA)

  # Optionally, filter for PP3 or continue with processing
  # df <- df %>% filter(ACMG_PP3 == "PP3")

  # Clean up, removing temporary columns used for calculation
  df <- df %>% select(-ACMG_PP3_count, -starts_with("cond_"))

  return(df)
}







# Process master function ----
# Master function to process genetic data
# Assuming you've already defined process_genetic_data and other required functions above
process_genetic_data <- function(input_path, samples_file_path, af_threshold) {
  imported_data <- read_data_file(input_path, samples_file_path)
  processed_data_list <- list()

  for (file_name in names(imported_data)) {
    cat("\nProcessing file: ", file_name, "\n")
    df <- imported_data[[file_name]]

    # Previous steps...
    df <- set_comp_het_flag(df)
    df <- filter_by_allele_frequency(df, af_threshold)
    df <- apply_acmg_pvs1(df)

    # Apply PS1 and PS5 criteria
    df <- apply_acmg_ps1(df)
    df <- apply_acmg_ps5(df)
    df <- apply_acmg_pm2(df)
    df <- apply_acmg_pm3(df)
    df <- apply_acmg_pp3(df)

    processed_data_list[[file_name]] <- df
  }

  return(processed_data_list)
}


# acmg tally  ----
# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(grid)
library(scico)

prepare_acmg_labels <- function(df) {
  acmg_labels <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5",
                   "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6",
                   "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4")

  # Ensure all ACMG columns exist in the dataframe, add them as NA if they do not
  for (acmg_label in acmg_labels) {
    if (!acmg_label %in% names(df)) {
      df[[acmg_label]] <- NA
    }
  }

  # Use coalesce to find the first non-NA ACMG label for each row, if needed
  # This step assumes you want to identify the highest priority ACMG classification per row
  # If you intended something different with ACMG_label, please adjust accordingly
  df$ACMG_highest <- coalesce(!!!select(df, all_of(acmg_labels)))

  # Count the number of non-NA ACMG labels for each row
  df$ACMG_count <- rowSums(!is.na(select(df, all_of(acmg_labels))), na.rm = TRUE)

  return(df)
}

# df <- prepare_acmg_labels(df)
# print(df)


# Plotting summaries ----
plot_criteria_count_each_gene <- function(df, file_suffix) {
  p <- df |>
    filter(ACMG_count > 1) |>
    ggplot(aes(y = ACMG_count, x = SYMBOL)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Gene symbol") +
    ylab("ACMG criteria count (>1)")

  ggsave(paste0("../output/", file_suffix, "_criteria_count_each_gene.pdf"), plot = p)
}

plot_criteria_gene_total <- function(df, file_suffix) {
  p <- df %>%
    group_by(SYMBOL) %>%
    summarise(acmg_count_per_symbol = sum(ACMG_count, na.rm = TRUE)) %>%
    na.omit() %>%
    ggplot(aes(x = acmg_count_per_symbol, fill = ..x..)) +
    geom_histogram(stat = "count", binwidth = 1, color = "black") +
    theme_minimal() +
    xlab("No. ACMG criteria (P) variants per gene") +
    ylab("Number of genes") +
    geom_text(stat = 'count', aes(label = ..count.., y = ..count.. + 20), color = "black") +
    guides(fill = FALSE) +
    scale_fill_scico(palette = 'acton', direction = 1)

  ggsave(paste0("../output/", file_suffix, "_criteria_gene_total.pdf"), plot = p)
}

plot_variants_per_criteria <- function(df, file_suffix) {
  p <- df |>
    ggplot(aes(x = ACMG_count, fill = ..x..)) +
    geom_histogram(binwidth = 1, color = "black") +
    xlab("No. ACMG criteria\nassigned (P)") +
    ylab("No. variants") +
    theme_minimal() +
    geom_text(stat = 'count', aes(label = ..count.., y = ..count.. + 300), color = "black") +
    guides(fill = FALSE) +
    scale_fill_scico(palette = 'acton', direction = 1)

  ggsave(paste0("../output/", file_suffix, "_variants_per_criteria.pdf"), plot = p, width = 9, height = 5)
}

# Sample usage
file_suffix <- "example_suffix"

# Assuming df is your dataframe
# df <- prepare_acmg_labels(df)

# plot_criteria_count_each_gene(df, file_suffix)
# plot_criteria_gene_total(df, file_suffix)
# plot_variants_per_criteria(df, file_suffix)

# Continue with other plotting functions similarly




#
#
#
#
#
#
