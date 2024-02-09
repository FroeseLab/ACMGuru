#' @title Apply ACMG PS5 Criterion
#'
#' @description Applies the ACMG PS5 criterion based on strong pathogenic evidence.
#'
#' @param df A dataframe containing genetic data.
#' @return A modified dataframe with the ACMG PS5 criterion applied.
#' @export
#' @import dplyr
apply_acmg_ps5 <- function(df) {
  cat("Applying ACMG PS5 criterion...\n")

  if ("IMPACT" %in% colnames(df)) {
    df <- df %>%
      dplyr::group_by(sample, SYMBOL) %>%
      dplyr::mutate(ACMG_PS5 = ifelse(any(IMPACT == "HIGH") & comp_het_flag == 1, "PS5", NA)) %>%
      dplyr::ungroup()
  } else {
    cat("Necessary columns for PS5 not found. Skipping PS5.\n")
  }

  return(df)
}

