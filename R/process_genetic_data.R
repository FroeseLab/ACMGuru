#' Process Genetic Data According to ACMG Guidelines
#'
#' This function processes genetic data by applying a series of ACMG criteria to identify variants based on specified thresholds and conditions.
#'
#' @param input_path Path to the input files or directory containing the CSV files.
#' @param samples_file_path Path to the samples file containing phenotype information.
#' @param af_threshold Allele frequency threshold for filtering variants.
#'
#' @return A list of data frames, each representing processed genetic data for a file.
#' @export
#' @importFrom dplyr filter
#' @importFrom tidyr separate
#' @importFrom stringr str_replace
#' @examples
#' input_path <- system.file("extdata", package = "YourPackageName")
#' samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")
#' af_threshold <- 0.01
#' processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)
process_genetic_data <- function(input_path, samples_file_path, af_threshold) {
  imported_data <- read_data_file(input_path, samples_file_path)
  processed_data_list <- list()

  for (file_name in names(imported_data)) {
    cat("\nProcessing file: ", file_name, "\n")
    df <- imported_data[[file_name]]


    # Apply ACMG criteria
    df <- set_comp_het_flag(df)
    df <- filter_by_allele_frequency(df, af_threshold)
    df <- apply_acmg_pvs1(df)
    df <- apply_acmg_ps1(df)
    df <- apply_acmg_ps5(df)
    df <- apply_acmg_pm2(df)
    df <- apply_acmg_pm3(df)
    df <- apply_acmg_pp3(df)

    processed_data_list[[file_name]] <- df
  }

  return(processed_data_list)
}

