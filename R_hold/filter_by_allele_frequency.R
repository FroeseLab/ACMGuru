#' Filter Data Based on Allele Frequency
#'
#' This function filters variants in the dataset based on a specified allele frequency threshold.
#'
#' @param df A dataframe containing genetic data, expected to have an `AF.x` column representing allele frequency.
#' @param af_threshold A numeric value specifying the maximum allele frequency for included variants.
#' @return A filtered dataframe with variants below the specified allele frequency threshold.
#' @import dplyr
#' @export
#' @examples
#' df <- data.frame(
#'   AF.x = c(0.01, 0.05, 0.2)
#' )
#' result <- filter_by_allele_frequency(df, 0.1)
filter_by_allele_frequency <- function(df, af_threshold) {
  df_filtered <- df %>%
    dplyr::filter(AF.x < af_threshold)
  
  return(df_filtered)
}

