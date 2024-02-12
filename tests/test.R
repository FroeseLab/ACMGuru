# Load necessary libraries
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(ggplot2)
# library(ggpubr)
# library(grid)
# library(data.table) # fread

# devtools::document()

# From the package root directory
library(devtools)
# devtools::install()
load_all() # Automatically points to the current directory as the package root

# varsome ----
# LE = less than equal to, GE = greater than equal to
# varsome <- read.csv(file = "../../data/singlecase/varsome_calibrated_insilico_thresholds.tsv", sep="\t")


# user variables ----
output_path <- "../output/"

samples_file_path <- "../data/samples.tsv" # phenotype data

input_path <- "../data/" # default: all files, or named file.
# input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"
# specific_files <- c("path/to/file1.csv", "path/to/file2.csv")
# input_path <- system.file("extdata", package = "YourPackageName")
# samples_file_path <- system.file("extdata", "samples.tsv", package = "YourPackageName")

file_list <- c(
  "../data/phrtmma_v1_chr21_40411318_41411317.csv",
  "../data/phrtmma_v1_chr21_41411318_42411317.csv",
  "../data/phrtmma_v1_chr21_42411318_43411317.csv"
)
processed_data_list <- process_genetic_data(NULL, samples_file_path, af_threshold, file_list)


af_threshold <- 0.1  # Allele frequency threshold
gnomad_max <- 1e-6
# metadataclass_file_path <- "../output/custom_reference_metadata.csv"

# start analysis ----
processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)



# process all files
input_path <- "../data/"
samples_file_path <- "../data/samples.tsv"
processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)

# single file
input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"
samples_file_path <- "../data/samples.tsv"
processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)

# list of files
input_path <- c(
  "../data/phrtmma_v1_chr21_40411318_41411317.csv",
  "../data/phrtmma_v1_chr21_41411318_42411317.csv",
  "../data/phrtmma_v1_chr21_42411318_43411317.csv"
)
samples_file_path <- "../data/samples.tsv"
processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)


# Check that all imported data chucks detected the correct column classes
# compare_column_classes_and_output_csv(processed_data_list, metadataclass_file_path)
compare_column_classes_and_output_csv(processed_data_list)

# Apply the corrected column classes to each dataset in processed_data_list
# processed_data_list <- apply_column_classes_to_processed_data(processed_data_list, metadataclass_file_path)

processed_data_list <- apply_column_classes_to_processed_data(processed_data_list)

# Merge all dataframes in the list into a single dataframe
all_data <- bind_rows(processed_data_list)

# processed_data_list[[1]] |> nrow()
# processed_data_list[[2]] |> nrow()
# processed_data_list[[1]] |> ncol()
# processed_data_list[[2]] |> ncol()

# df1 <- processed_data_list[[1]]
# df2 <- processed_data_list[[2]]

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
# write.csv(all_data, paste0("../output/aggregated_data_", Sys.Date(), ".csv"), row.names = FALSE)

