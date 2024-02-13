#' @title Apply ACMG PM2 Criterion
#'
#' @description Applies the ACMG PM2 criterion based on allele frequency in population databases.
#'
#' @param df A dataframe containing genetic data.
#' @param gnomad_max The maximum allele frequency considered for PM2.
#' @return A modified dataframe with the ACMG PM2 criterion applied.
#' @export
#' @import dplyr
apply_acmg_pm2 <- function(df, gnomad_max = 1e-6) {
  cat("Applying ACMG PM2 criterion...\n")

  if ("gnomAD_AF" %in% colnames(df)) {
    df <- df %>%
      dplyr::mutate(ACMG_PM2 = ifelse(gnomAD_AF < gnomad_max, "PM2", NA))
  } else {
    cat("Column 'gnomAD_AF' not found. Skipping PM2.\n")
  }

  return(df)
}

