#' Plotting Summaries for ACMG Criteria
#'
#' Contains functions for generating various plots to summarize the counts and distributions of ACMG criteria across genes and variants.

#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import scico

# Plot Criteria Count for Each Gene
plot_criteria_count_each_gene <- function(df, file_suffix) {
  p <- df %>%
    dplyr::filter(ACMG_count > 1) %>%
    ggplot2::ggplot(ggplot2::aes(y = ACMG_count, x = SYMBOL)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Gene symbol") +
    ggplot2::ylab("ACMG criteria count (>1)")

  ggplot2::ggsave(paste0("../output/", file_suffix, "_criteria_count_each_gene.pdf"), plot = p)
}

# Plot Criteria Gene Total
plot_criteria_gene_total <- function(df, file_suffix) {
  p <- df %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(acmg_count_per_symbol = sum(ACMG_count, na.rm = TRUE)) %>%
    na.omit() %>%
    ggplot2::ggplot(ggplot2::aes(x = acmg_count_per_symbol, fill = ..x..)) +
    ggplot2::geom_histogram(stat = "count", binwidth = 1, color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::xlab("No. ACMG criteria (P) variants per gene") +
    ggplot2::ylab("Number of genes") +
    ggplot2::geom_text(stat = 'count', ggplot2::aes(label = ..count.., y = ..count.. + 20), color = "black") +
    ggplot2::guides(fill = FALSE) +
    scico::scale_fill_scico(palette = 'acton', direction = 1)

  ggplot2::ggsave(paste0("../output/", file_suffix, "_criteria_gene_total.pdf"), plot = p)
}

# Plot Variants Per Criteria
plot_variants_per_criteria <- function(df, file_suffix) {
  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = ACMG_count, fill = ..x..)) +
    ggplot2::geom_histogram(binwidth = 1, color = "black") +
    ggplot2::xlab("No. ACMG criteria\nassigned (P)") +
    ggplot2::ylab("No. variants") +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(stat = 'count', ggplot2::aes(label = ..count.., y = ..count.. + 300), color = "black") +
    ggplot2::guides(fill = FALSE) +
    scico::scale_fill_scico(palette = 'acton', direction = 1)

  ggplot2::ggsave(paste0("../output/", file_suffix, "_variants_per_criteria.pdf"), plot = p, width = 9, height = 5)
}

#' @examples
#' file_suffix <- "example_suffix"
#'
#' # Assuming df is your dataframe prepared with ACMG labels
#' # df <- prepare_acmg_labels(df)
#'
#' # plot_criteria_count_each_gene(df, file_suffix)
#' # plot_criteria_gene_total(df, file_suffix)
#' # plot_variants_per_criteria(df, file_suffix)

