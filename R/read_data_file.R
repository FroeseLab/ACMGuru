#' Read Data Files and Merge with Sample Phenotype Information
#'
#' This function reads data files from a specified path and merges them with sample phenotype information.
#'
#' @param input_path Path to the input files or directory containing the files.
#' @param samples_file_path Path to the samples file containing phenotype information.
#'
#' @return A list of data frames, each representing a merged data file.
#' @import data.table
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#' input_path <- system.file("extdata", package = "YourPackageName")
#' samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")
#' read_data_file(input_path, samples_file_path)
read_data_file <- function(input_path, samples_file_path) {
  files <- character()

  if (dir.exists(input_path)) {
    files <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
  } else if (file.exists(input_path) && grepl("\\.csv$", input_path)) {
    files <- input_path
  } else {
    stop("Input path is neither a valid file nor a directory.")
  }

  cat("Reading samples information from: ", samples_file_path, "\n")
  samples_df <- fread(samples_file_path, sep = "\t", col.names = c("sample", "cohort_pheno"), header = FALSE, data.table = FALSE)

  data_list <- list()
  for (file_path in files) {
    cat("Reading file: ", file_path, "\n")
    df <- fread(file_path, sep = ",", data.table = FALSE) # fread reads into a data.table by default, data.table = FALSE converts it to a data.frame
    df <- merge(df, samples_df, by = "sample", all.x = TRUE)
    data_list[[basename(file_path)]] <- df
  }

  return(data_list)
}


