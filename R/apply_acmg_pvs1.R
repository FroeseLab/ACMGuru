#' @title Apply ACMG PVS1 Criterion
#'
#' @description This function applies the ACMG PVS1 criterion to a given dataframe.
#'
#' @param df A dataframe containing genetic data.
#' @return A modified dataframe with the ACMG PVS1 criterion applied.
#' @export
#' @import dplyr
apply_acmg_pvs1 <- function(df) {
  cat("Applying ACMG PVS1 criterion...\n")
  
  df <- df %>%
    dplyr::mutate(ACMG_PVS1 = NA) %>%
    dplyr::mutate(ACMG_PVS1 = ifelse(IMPACT == "HIGH" & genotype == 2, "PVS1", ACMG_PVS1))
  
  if ("Inheritance" %in% colnames(df)) {
    df <- df %>%
      dplyr::mutate(ACMG_PVS1 = ifelse(IMPACT == "HIGH" & Inheritance == "AD", "PVS1", ACMG_PVS1))
  }
  
  df <- df %>%
    dplyr::group_by(sample, SYMBOL) %>%
    dplyr::mutate(ACMG_PVS1 = ifelse(ACMG_PVS1 == "PVS1", "PVS1",
                                     ifelse(sum(IMPACT == "HIGH" & comp_het_flag == 1, na.rm = TRUE) >= 2 & IMPACT == "HIGH", "PVS1", ACMG_PVS1))) %>%
    dplyr::ungroup()
  
  return(df)
}

