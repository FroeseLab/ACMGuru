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
    
    processed_data_list[[file_name]] <- df
  }
  
  return(processed_data_list)
}
