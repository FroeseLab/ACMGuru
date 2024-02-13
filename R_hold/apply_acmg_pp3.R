#' Apply ACMG PP3 Criterion
#'
#' Applies the ACMG PP3 criterion based on multiple lines of computational evidence supporting a deleterious effect on the gene or gene product. This function assesses variants for evidence of pathogenicity based on computational predictions.
#'
#' @param df A dataframe containing genetic data with necessary predictive columns.
#' @return A modified dataframe with the ACMG PP3 criterion applied.
#' @export
#' @import dplyr
#' @import tidyr
#' @import stringr
apply_acmg_pp3 <- function(df) {
  cat("Applying ACMG PP3 criterion...\n")

  # Check if the necessary columns exist
  required_columns <- c("CADD_PHRED", "REVEL_rankscore", "MetaLR_pred",
                        "MutationAssessor_pred", "SIFT", "PolyPhen")
  if(!all(required_columns %in% colnames(df))) {
    cat("One or more required columns for PP3 are missing. Skipping PP3.\n")
    return(df)
  }

  # PP3 ----
  # Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.). CUSTOM: assigned if >=3 thresholds passed.
  # In-Silico Predictions: VarSome now implements the ClinGen recommendations from Evidence-based calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for clinical use of PP3/BP4 criteria: # Only one engine at a time is used, depending on availability of data, in order: MitoTip & MitImpact, MetaRNN, CADD (Premium only), DANN (if CADD is not available). # The maximum strength allowed for rules PP3 & BP4 is Strong, even if there may be evidence for Very Strong, with the exception of variants that are predicted splicing (ie: similar to PVS1). # The strength is limited to Supporting, if there's Moderate evidence from rules PM1 or PM5. # Splice prediction (scSNV) is given priority over the other in-silico predictions. # conservation is used for some low-sensitivity variant types, or if no other in-silico prediction is available. Please refer to PP3 and BP4 for more specific detail.

  # Prepare data for conditions
  df <- df %>%
    tidyr::separate(SIFT, into = c("SIFT_label", "SIFT_score"), sep = "\\(", remove = TRUE) %>%
    dplyr::mutate(SIFT_score = stringr::str_replace(SIFT_score, "\\)", "")) %>%
    tidyr::separate(PolyPhen, into = c("PolyPhen_label", "PolyPhen_score"), sep = "\\(", remove = TRUE) %>%
    dplyr::mutate(PolyPhen_score = stringr::str_replace(PolyPhen_score, "\\)", "")) %>%
    dplyr::mutate(
      CADD_PHRED = as.numeric(CADD_PHRED),
      REVEL_rankscore = as.numeric(REVEL_rankscore),
      ACMG_PP3 = NA
    )

  # Define conditions and count them
  df$ACMG_PP3_count <- rowSums(cbind(
    df$CADD_PHRED >= 30,
    df$REVEL_rankscore > 0.5,
    df$MetaLR_pred == "D",
    df$MutationAssessor_pred == "H",
    df$SIFT_label == "deleterious",
    df$PolyPhen_label == "probably_damaging"
  ), na.rm = TRUE)

  # Apply PP3 if the threshold is met
  df$ACMG_PP3 <- ifelse(df$ACMG_PP3_count >= 3, "PP3", NA)

  # Optionally, filter for PP3 or continue with processing
  # df <- df %>% dplyr::filter(ACMG_PP3 == "PP3")

  # Clean up, removing temporary columns used for calculation
  df <- df %>% dplyr::select(-ACMG_PP3_count, -starts_with("cond_"))

  return(df)
}

