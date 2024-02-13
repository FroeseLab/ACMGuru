#' Apply Column Classes Based on Reference Metadata
#'
#' Adjusts the column classes of a list of data frames based on a provided reference metadata CSV file.
#'
#' @param processed_data_list A list of data frames to adjust.
#' @param metadataclass_file_path Path to the CSV file containing reference metadata for column classes.
#'
#' @return A list of adjusted data frames with updated column classes.
#' @import utils
#' @export
#' @examples
#' processed_data_list <- list(df1 = data.frame(a = as.character(1:3)), df2 = data.frame(a = as.character(4:6)))
#' metadataclass_file_path <- system.file("extdata", "metadata.csv", package = "YourPackageName")
#' apply_column_classes(processed_data_list, metadataclass_file_path)
apply_column_classes <- function(processed_data_list, metadataclass_file_path) {
  if (!file.exists(metadataclass_file_path)) {
    cat("Reference metadata file does not exist. Skipping column class adjustment.\n")
    return(processed_data_list)
  }

  reference_metadata <- read.csv(metadataclass_file_path, stringsAsFactors = FALSE)
  reference_metadata_list <- setNames(as.list(reference_metadata$most_common_class), reference_metadata$column_name)

  adjusted_data_list <- lapply(processed_data_list, function(df) {
    for (column_name in names(reference_metadata_list)) {
      if (column_name %in% names(df)) {
        class_type <- reference_metadata_list[[column_name]]
        df[[column_name]] <- tryCatch({
          if (class_type == "factor") as.factor(df[[column_name]])
          else if (class_type == "numeric") as.numeric(df[[column_name]])
          else if (class_type == "integer") as.integer(df[[column_name]])
          else if (class_type == "character") as.character(df[[column_name]])
          else df[[column_name]] # Return column as-is if class type not recognized
        }, error = function(e) {
          cat("Error converting column", column_name, "to", class_type, "\n")
          df[[column_name]] # Return column as-is in case of error
        })
      }
    }
    return(df)
  })

  return(adjusted_data_list)
}

