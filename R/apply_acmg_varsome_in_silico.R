#' @title Apply Varsome In-silico Predictions for ACMG PP3 Criterion
#'
#' @description Enhances the ACMG PP3 criterion by incorporating Varsome in-silico
#' prediction scores to assess variant pathogenicity.
#'
#' @param df A dataframe containing genetic data with necessary predictive columns.
#' @param varsome_data A dataframe containing Varsome in-silico prediction scores and thresholds.
#' @param pathogenic_types A character vector of pathogenic types to assess ("Strong_pathogenic_GE",
#' "Moderate_pathogenic_GE", "Supporting_pathogenic_GE").
#' @return A dataframe with additional columns indicating the count of pathogenic predictions met for each type.
#' @export
#' @import dplyr
#' @import tidyr
#' @import stringr
# apply_varsome_in_silico <- function(df, varsome_file_path = NULL, pathogenic_types = c("Strong_pathogenic_GE", "Moderate_pathogenic_GE", "Supporting_pathogenic_GE")) {
#   # Use the default file if no path is provided
#   if (is.null(varsome_file_path)) {
#     varsome_file_path <- system.file("extdata", "varsome_calibrated_insilico_thresholds.tsv", package = "ACMGuru")
#   }
#
#   # Check if the file exists
#   if (!file.exists(varsome_file_path)) {
#     stop("Varsome data file does not exist.")
#   }
#
#   # Load varsome_data
#   varsome_data <- fread(varsome_file_path, sep = "\t")
#
#   # Ensure varsome_data is properly formatted and contains necessary columns
#   required_columns <- c("Engine", "Strong_pathogenic_GE", "Moderate_pathogenic_GE", "Supporting_pathogenic_GE")
#   if (!all(required_columns %in% colnames(varsome_data))) {
#     stop("Varsome data is missing required columns.")
#   }
#
#   # Prepare Varsome data: match engine names with df columns
#   # Define replacements as a named vector
#   names_to_replace <- c(
#     BayesDel_addAF = "BayesDel_addAF_score",
#     BayesDel_noAF = "BayesDel_noAF_score",
#     CADD = "CADD_PHRED",
#     DANN = "DANN_score",
#     EIGEN = "Eigen.raw_coding",
#     `EIGEN-PC` = "Eigen.PC.phred_coding",
#     FATHMM = "FATHMM_score",
#     `FATHMM-MKL` = "fathmm.MKL_coding_score",
#     `FATHMM-XF` = "fathmm.XF_coding_score",
#     LRT = "LRT_score",
#     `M-CAP` = "M.CAP_score",
#     MetaLR = "MetaLR_score",
#     MetaSVM = "MetaSVM_score",
#     MetaRNN = "MetaRNN_score",
#     MutPred = "MutPred_score",
#     MutationAssessor = "MutationAssessor_score",
#     MutationTaster = "MutationTaster_score",
#     phastCons100way_vertebrate = "phastCons100way_vertebrate",
#     `Polyphen2-HDIV` = "Polyphen2_HDIV_score",
#     `Polyphen2-HVAR` = "Polyphen2_HVAR_score",
#     PROVEAN = "PROVEAN_score",
#     REVEL = "REVEL_score",
#     SIFT = "SIFT_score"
#   )
#
#   # Apply the renaming
#   for (old_name in names(names_to_replace)) {
#     new_name <- names_to_replace[old_name]
#     if (old_name %in% colnames(varsome_data)) {
#       varsome_data <- rename(varsome_data, !!new_name := !!sym(old_name))
#     }
#   }
#
#
#   # Apply the scores based on Varsome thresholds
#   for (pathogenic_type in pathogenic_types) {
#     df <- calculate_varsome_score(df, varsome_data, pathogenic_type)
#   }
#
#   return(df)
# }

#' Calculate Varsome Score
#'
#' @description Internal function to calculate scores based on Varsome in-silico predictions.
# calculate_varsome_score <- function(df, varsome_data, pathogenic_type) {
#   varsome_list <- varsome_data %>%
#     filter(pathogenicity == pathogenic_type) %>%
#     select(Engine, threshold)
#
#   df[[pathogenic_type]] <- 0
#
#   for (i in 1:nrow(varsome_list)) {
#     engine <- varsome_list$Engine[i]
#     threshold <- varsome_list$threshold[i]
#
#     if (!(engine %in% names(df))) {
#       next # Skip if the engine is not found in the dataframe
#     }
#
#     # Convert to numeric if not already
#     df[[engine]] <- as.numeric(df[[engine]])
#
#     # Add to pathogenic_type count if the score meets or exceeds the threshold
#     df[[pathogenic_type]] <- df[[pathogenic_type]] + (df[[engine]] >= threshold)
#   }
#
#   return(df)
# }

Strong_pathogenic_GE_threshold <- 3
Moderate_pathogenic_GE_threshold <- 4
Supporting_pathogenic_GE_threshold <- 10

apply_varsome_in_silico <- function(df, varsome_file_path = NULL, pathogenic_types = c("Strong_pathogenic_GE", "Moderate_pathogenic_GE", "Supporting_pathogenic_GE")) {
  # Default file path if not provided
  if (is.null(varsome_file_path)) {
    varsome_file_path <- system.file("extdata", "varsome_calibrated_insilico_thresholds.tsv", package = "ACMGuru")
  }

  # Load varsome_data
  varsome_data <- fread(varsome_file_path, sep = "\t")

  # Ensure varsome_data has required columns
  if (!all(c("Engine", "Strong_pathogenic_GE", "Moderate_pathogenic_GE", "Supporting_pathogenic_GE") %in% colnames(varsome_data))) {
    stop("Varsome data is missing required columns.")
  }

  # Renaming logic: Match engine names with df columns
  names_to_replace <- list(
    c("BayesDel_addAF", "BayesDel_addAF_score"),
    c("BayesDel_noAF", "BayesDel_noAF_score"),
    c("CADD", "CADD_PHRED"),
    c("DANN", "DANN_score"),
    c("EIGEN", "Eigen.raw_coding"),
    c("EIGEN-PC", "Eigen.PC.phred_coding"),
    c("FATHMM", "FATHMM_score"),
    c("FATHMM-MKL", "fathmm.MKL_coding_score"),
    c("FATHMM-XF", "fathmm.XF_coding_score"),
    c("LRT", "LRT_score"),
    c("M-CAP", "M.CAP_score"),
    c("MetaLR", "MetaLR_score"),
    c("MetaSVM", "MetaSVM_score"),
    c("MetaRNN", "MetaRNN_score"),
    c("MutPred", "MutPred_score"),
    c("MutationAssessor", "MutationAssessor_score"),
    c("MutationTaster", "MutationTaster_score"),
    c("phastCons100way_vertebrate", "phastCons100way_vertebrate"),
    c("Polyphen2-HDIV", "Polyphen2_HDIV_score"),
    c("Polyphen2-HVAR", "Polyphen2_HVAR_score"),
    c("PROVEAN", "PROVEAN_score"),
    c("REVEL", "REVEL_score"),
    c("SIFT", "SIFT_score")
  )

  # varsome_data <- varsome_data %>%
  #   rename_at(vars(names(names_to_replace)), ~ names_to_replace[.])
  for (name_pair in names_to_replace) {
    varsome_data$Engine <- replace(varsome_data$Engine, varsome_data$Engine == name_pair[1], name_pair[2])
  }

  # Apply scores based on Varsome thresholds
  for (pathogenic_type in pathogenic_types) {
    df <- calculate_varsome_score(df, varsome_data, pathogenic_type)
  }

  return(df)
}

# calculate_varsome_score <- function(df, varsome_data, pathogenic_type) {
#   # Initialize the count column for the pathogenic type
#   df[[pathogenic_type]] <- 0
#
#   # Iterate through each engine in varsome_data
#   engines <- unique(varsome_data$Engine)
#   for (engine in engines) {
#     threshold <- varsome_data[varsome_data$Engine == engine & varsome_data$pathogenicity == pathogenic_type, "threshold"]
#
#     # Apply threshold if engine exists in df
#     if (engine %in% names(df)) {
#       condition <- df[[engine]] >= threshold
#       df[[pathogenic_type]] <- df[[pathogenic_type]] + as.integer(condition)
#     }
#   }
#
#   return(df)
# }

calculate_varsome_score <- function(df, varsome_data, pathogenic_type) {
  # Initialize the count column for the pathogenic type
  df[[pathogenic_type]] <- 0

  # Ensure varsome_data includes thresholds for the specified pathogenic type
  if (!pathogenic_type %in% names(varsome_data)) {
    stop(paste("Pathogenic type", pathogenic_type, "not found in varsome_data columns."))
  }

  # Iterate through each engine in varsome_data
  engines <- unique(varsome_data$Engine)
  for (engine in engines) {
    # Directly access the threshold for the engine and pathogenic type
    # Assuming varsome_data is structured with one row per engine, and columns for each pathogenic type's threshold
    engine_row <- varsome_data[varsome_data$Engine == engine, ]
    if (nrow(engine_row) == 0) {
      next  # Skip if engine is not found in varsome_data
    }

    # Extract the threshold for the current pathogenic type
    threshold <- engine_row[[pathogenic_type]]

    # Apply threshold if engine exists in df
    if (engine %in% names(df)) {
      condition <- df[[engine]] >= threshold
      df[[pathogenic_type]] <- df[[pathogenic_type]] + as.integer(condition)
    }
  }
  cat(paste0("Applying threshold: ",pathogenic_type,"\n"))
  return(df)
}
