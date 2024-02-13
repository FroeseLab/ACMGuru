#' Read Data Files and Merge with Sample Phenotype Information
#'
#' This function reads data files from a specified path and merges them with sample phenotype information. It can handle a single file, all files within a directory, or a specific list of files.
#'
#' @param input_path Optional. Path to the input directory containing the files or a single file path. Used only if file_list is NULL.
#' @param samples_file_path Path to the samples file containing phenotype information.
#' @param file_list Optional. A vector of specific file paths to be processed. Overrides input_path if provided.
#' @return A list of data frames, each representing a merged data file.
#' @importFrom data.table fread
#' @export
#' @examples
# Process all files in a directory:
# input_path <- system.file("extdata", package = "YourPackageName")
# samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")
# processed_data_list <- read_data_file(input_path, samples_file_path)

# To process a single file:
# single_file <- system.file("extdata", "single_file.csv", package = "YourPackageName")
# samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")
# processed_data_list <- read_data_file(single_file, samples_file_path)

# To process a specific list of files:
# specific_files <- c(
#   system.file("extdata", "file1.csv", package = "YourPackageName"),
#   system.file("extdata", "file2.csv", package = "YourPackageName")
# )
# samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")
# processed_data_list <- read_data_file(NULL, samples_file_path, specific_files)

read_data_file <- function(input_path = NULL, samples_file_path, file_list = NULL) {
  files <- character()

  # Check if a specific list of files is provided
  if (!is.null(file_list) && length(file_list) > 0) {
    files <- file_list
  } else if (!is.null(input_path)) {
    if (length(input_path) == 1) {
      # Single input_path provided
      if (dir.exists(input_path)) {
        # It's a directory
        files <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
      } else if (file.exists(input_path) && grepl("\\.csv$", input_path)) {
        # It's a single file
        files <- input_path
      } else {
        stop("Input path is neither a valid file nor a directory.")
      }
    } else if (length(input_path) > 1) {
      # input_path is treated as a list of files
      files <- input_path
    } else {
      stop("Please provide either a valid directory path, a single file path, or a list of file paths.")
    }
  } else {
    stop("input_path or file_list must be provided.")
  }

  cat("Reading samples information from: ", samples_file_path, "\n")
  samples_df <- fread(samples_file_path, sep = "\t", col.names = c("sample", "cohort_pheno"), header = FALSE, data.table = FALSE)

  data_list <- list()
  for (file_path in files) {
    cat("Reading file: ", file_path, "\n")
    df <- fread(file_path, sep = ",", data.table = FALSE)
    df <- merge(df, samples_df, by = "sample", all.x = TRUE)
    data_list[[basename(file_path)]] <- df
  }

  return(data_list)
}
