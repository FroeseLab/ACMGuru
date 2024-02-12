#' Read Data Files and Merge with Sample Phenotype Information
#'
#' @param input_path Optional. Path to the input directory containing the files or a single file path.
#'                   Used only if file_list is NULL.
#' @param samples_file_path Path to the samples file containing phenotype information.
#' @param file_list Optional. A vector of specific file paths to be processed. Overrides input_path if provided.
#' @return A list of data frames, each representing a merged data file.
#' @importFrom data.table fread
#' @export
#' @examples
#' input_path <- system.file("extdata", package = "YourPackageName")
#' samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")
#' specific_files <- c("path/to/file1.csv", "path/to/file2.csv")
#' read_data_file(input_path, samples_file_path, specific_files)
read_data_file <- function(input_path = NULL, samples_file_path, file_list = NULL) {
  files <- character()

  if (!is.null(file_list) && length(file_list) > 0) {
    # Use the specific list of files provided
    files <- file_list
  } else if (!is.null(input_path) && dir.exists(input_path)) {
    # If input_path is a directory and file_list is not provided
    files <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
  } else if (!is.null(input_path) && file.exists(input_path) && grepl("\\.csv$", input_path)) {
    # If input_path is a single file and file_list is not provided
    files <- input_path
  } else {
    stop("Please provide either a valid directory path, a single file path, or a list of file paths.")
  }

  cat("Reading samples information from: ", samples_file_path, "\n")
  samples_df <- fread(samples_file_path, sep = "\t", col.names = c("sample", "cohort_pheno"), header = FALSE, data.table = FALSE)

  data_list <- list()
  for (file_path in files) {
    cat("Reading file: ", file_path, "\n")
    df <- fread(file_path, sep = ",", data.table = FALSE)  # fread reads into a data.table by default, data.table = FALSE converts it to a data.frame
    df <- merge(df, samples_df, by = "sample", all.x = TRUE)
    data_list[[basename(file_path)]] <- df
  }

  return(data_list)
}
