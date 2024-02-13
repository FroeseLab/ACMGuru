#' @title Apply ACMG PS1 Criterion
#'
#' @description This function applies the ACMG PS1 criterion, identifying pathogenic variants based on clinical significance.
#'
#' @param df A dataframe containing genetic data with a `CLIN_SIG` column.
#' @return A modified dataframe with the ACMG PS1 criterion applied.
#' @export
#' @import dplyr
apply_acmg_ps1 <- function(df) {
  cat("Applying ACMG PS1 criterion...\n")
  
  if ("CLIN_SIG" %in% colnames(df)) {
    df <- df %>%
      dplyr::mutate(ACMG_PS1 = ifelse(CLIN_SIG == "pathogenic", "PS1", NA))
  } else {
    cat("Column 'CLIN_SIG' not found. Skipping PS1.\n")
  }
  
  return(df)
}

