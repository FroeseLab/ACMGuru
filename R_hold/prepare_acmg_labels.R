#' Prepare ACMG Labels
#'
#' Prepares ACMG labels by ensuring all expected ACMG label columns exist in the dataframe, assigning NA where they do not, identifying the highest priority ACMG classification for each row, and counting the number of non-NA ACMG labels for each row.
#'
#' @param df A dataframe containing genetic data potentially including various ACMG label columns.
#' @return A modified dataframe with additional columns for the highest priority ACMG classification (`ACMG_highest`) and the count of non-NA ACMG labels (`ACMG_count`) for each row.
#' @export
#' @import dplyr
#' @import tidyr
#' @examples
#' df <- data.frame(
#'   ACMG_PVS1 = c(NA, "PVS1"),
#'   ACMG_PS1 = c("PS1", NA)
#' )
#' df_prepared <- prepare_acmg_labels(df)
prepare_acmg_labels <- function(df) {
  acmg_labels <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5",
                   "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6",
                   "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4")

  # Ensure all ACMG columns exist in the dataframe, add them as NA if they do not
  for (acmg_label in acmg_labels) {
    if (!acmg_label %in% names(df)) {
      df[[acmg_label]] <- NA
    }
  }

  # Use dplyr::coalesce to find the first non-NA ACMG label for each row
  df$ACMG_highest <- dplyr::coalesce(!!!dplyr::select(df, dplyr::all_of(acmg_labels)))

  # Count the number of non-NA ACMG labels for each row
  df$ACMG_count <- rowSums(!is.na(dplyr::select(df, dplyr::all_of(acmg_labels))), na.rm = TRUE)

  return(df)
}

