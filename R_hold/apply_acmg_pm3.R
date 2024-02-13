#' @title Apply ACMG PM3 Criterion
#'
#' @description Applies the ACMG PM3 criterion for recessive disorders detected in trans with a pathogenic variant.
#'
#' @param df A dataframe containing genetic data.
#' @return A modified dataframe with the ACMG PM3 criterion applied.
#' @export
#' @import dplyr
apply_acmg_pm3 <- function(df) {
  cat("Applying ACMG PM3 criterion...\n")

  if (all(c("comp_het_flag", "ACMG_PS1", "ACMG_PS5") %in% colnames(df))) {
    df <- df %>%
      dplyr::group_by(sample, SYMBOL) %>%
      dplyr::mutate(ACMG_PM3 = ifelse(comp_het_flag == 1 & (ACMG_PS1 == "PS1" | ACMG_PS5 == "PS5"), "PM3", NA)) %>%
      dplyr::ungroup()
  } else {
    cat("Necessary columns for PM3 not found. Skipping PM3.\n")
  }

  return(df)
}

