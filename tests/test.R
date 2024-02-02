# varsome ----
# LE = less than equal to, GE = greater than equal to
varsome <- read.csv(file = "../../data/singlecase/varsome_calibrated_insilico_thresholds.tsv", sep="\t")


source("../R/acmguru_vcurrent.R")

if (!dir.exists("../output/")) {
  dir.create("../output/", recursive = TRUE)
}

# Specify the path for input data and samples file
input_path <- "../data/"
# input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"
samples_file_path <- "../data/samples.tsv"
af_threshold <- 0.1  # Allele frequency threshold
gnomad_max <- 1e-6

# Assuming process_genetic_data function returns a list of processed dataframes
processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)

# Merge all dataframes in the list into a single dataframe
# This assumes that all dataframes have the same structure

all_data <- bind_rows(processed_data_list)
processed_data_list[[1]] |> nrow()
processed_data_list[[2]] |> nrow()
processed_data_list[[1]] |> ncol()
processed_data_list[[2]] |> ncol()

processed_data_list[[1]]$BayesDel_addAF_pred
processed_data_list[[4]]$BayesDel_addAF_pred

# debug -----
# Assuming processed_data_list is your list of data frames

# # Function to get column types for a data frame
# get_column_types <- function(df) {
#   sapply(df, class)
# }
#
# # Initialize a list to hold column types for each data frame
# column_types_list <- lapply(processed_data_list, get_column_types)
#
# # Function to compare the column types across all data frames
# compare_column_types <- function(column_types_list) {
#   # Use the first item as the base for comparison
#   base_types <- column_types_list[[1]]
#   differing_columns <- list()
#
#   # Iterate over the list, comparing each set of types to the base
#   for (i in 2:length(column_types_list)) {
#     current_types <- column_types_list[[i]]
#     diff <- !(base_types %in% current_types) | !(current_types %in% base_types)
#     if (any(diff)) {
#       differing_columns[[length(differing_columns) + 1]] <- names(current_types[diff])
#     }
#   }
#   return(unique(unlist(differing_columns)))
# }
#
# # Get the names of columns with differing types
# columns_with_differing_types <- compare_column_types(column_types_list)
# print(columns_with_differing_types)
#
# # Optionally, print out the specific types for a column across data frames for investigation
# column_to_inspect <- "BayesDel_addAF_pred" # Example column name
# if(column_to_inspect %in% columns_with_differing_types) {
#   types_per_df <- sapply(column_types_list, function(types) types[column_to_inspect])
#   print(types_per_df)
# }

# debug ----

# Ensure all_data is indeed a dataframe
if (!is.data.frame(all_data)) {
  stop("all_data is not a dataframe")
}

# Now, you have a single dataframe 'all_data' with data from all processed files

# Apply ACMG labels preparation on the aggregated data
all_data <- prepare_acmg_labels(all_data)

# Define a generic suffix or use a specific identifier for your output files
file_suffix <- "aggregate"

# Generate plots based on the aggregated data
plot_criteria_count_each_gene(all_data, file_suffix)
plot_criteria_gene_total(all_data, file_suffix)
plot_variants_per_criteria(all_data, file_suffix)

# Optionally, you can save the aggregated dataframe for further analysis
write.csv(all_data, paste0("../output/aggregated_data_", Sys.Date(), ".csv"), row.names = FALSE)


# source("../R/acmguru_vcurrent.R")
#
# if (!dir.exists("../output/")) {
#   dir.create("../output/", recursive = TRUE)
# }
#
# # Specify the path for input data and samples file
# input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"
# samples_file_path <- "../data/samples.tsv"
# af_threshold <- 0.1  # Allele frequency threshold
#
# # Process the genetic data
# processed_data <- process_genetic_data(input_path, samples_file_path, af_threshold)
#
# # Example to view processed data of a specific file
# if (length(processed_data) > 0) {
#   specific_file_name <- names(processed_data)[1]  # Use the first file name as key
#   df_test <- processed_data[[specific_file_name]]  # Correctly access the dataframe
#
#   # Prepare ACMG labels for the test dataframe
#   df_test <- prepare_acmg_labels(df_test)
#
#   # Define a file suffix for naming the output files, extracted from the specific file name
#   file_suffix <- gsub(".csv$", "", basename(specific_file_name))
#
#   # Generate plots for this specific dataframe
#   plot_criteria_count_each_gene(df_test, file_suffix)
#   plot_criteria_gene_total(df_test, file_suffix)
#   plot_variants_per_criteria(df_test, file_suffix)
#
#   # Optionally, print the head of the processed dataframe to check the results
#   # print(head(df_test))
# }
