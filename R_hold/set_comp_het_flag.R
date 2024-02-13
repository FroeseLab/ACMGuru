#' Set Compound Heterozygous Flags
#'
#' This function sets flags for compound heterozygous variants within the dataset.
#' WARNING: This function does not check phase.
#'
#' @param df A dataframe containing genetic data with columns `sample`, `SYMBOL`, and `genotype`.
#' @return A modified dataframe with a new `comp_het_flag` column.
#' @import dplyr
#' @export
#' @examples
#' df <- data.frame(
#'   sample = c("Sample1", "Sample1"),
#'   SYMBOL = c("Gene1", "Gene1"),
#'   genotype = c(1, 2)
#' )
#' result <- set_comp_het_flag(df)
set_comp_het_flag <- function(df) {
  cat("Setting compound heterozygous flag...\n")
  df <- df %>%
    dplyr::group_by(sample, SYMBOL) %>%
    dplyr::mutate(comp_het_flag = ifelse(n() > 1, 1, NA)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(comp_het_flag = ifelse(is.na(comp_het_flag) & genotype == 2, 1, comp_het_flag))
  
  return(df)
}

