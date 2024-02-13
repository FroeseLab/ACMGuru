#' Apply Column Classes to Processed Data
#'
#' Applies the column classes to each dataframe in a list based on a provided reference metadata CSV file.
#'
#' @param processed_data_list A list of dataframes to be adjusted.
#' @param metadataclass_file_path Optional. Path to the CSV file containing reference metadata for column classes.
#' If not provided, uses the default file included in the package.
#' @export
apply_column_classes_to_processed_data <- function(processed_data_list, metadataclass_file_path = NULL) {
  # Use the default file if no path is provided
  if (is.null(metadataclass_file_path)) {
    metadataclass_file_path <- system.file("extdata", "suggested_reference_metadata.csv", package = "ACMGuru")
  }

  # Check if the file exists
  if (!file.exists(metadataclass_file_path)) {
    stop("Reference metadata file does not exist.")
  }

  reference_metadata <- read.csv(metadataclass_file_path, stringsAsFactors = FALSE)
  reference_metadata_list <- setNames(as.list(reference_metadata$most_common_class), reference_metadata$column_name)

  adjusted_data_list <- lapply(processed_data_list, function(df) {
    for (column_name in names(reference_metadata_list)) {
      if (column_name %in% names(df)) {
        class_type <- reference_metadata_list[[column_name]]
        df[[column_name]] <- tryCatch({
          switch(class_type,
                 "factor" = as.factor(df[[column_name]]),
                 "numeric" = as.numeric(df[[column_name]]),
                 "integer" = as.integer(df[[column_name]]),
                 "character" = as.character(df[[column_name]]),
                 df[[column_name]])  # Default case
        }, error = function(e) {
          message(paste("Error converting column", column_name, "to", class_type, "in dataframe. Error:", e$message))
          df[[column_name]]
        })
      }
    }
    return(df)
  })

  return(adjusted_data_list)
}
