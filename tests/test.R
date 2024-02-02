source("../R/acmguru_vcurrent.R")

# input_path <- "../data"
input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"
samples_file_path <- "../data/samples.tsv"
af_threshold <- 0.1  # Define the allele frequency threshold

# Process the genetic data
processed_data <- process_genetic_data(input_path, samples_file_path, af_threshold)

# Example to view processed data of a specific file
if (length(processed_data) > 0) {
    specific_file_name <- names(processed_data)[1]  # Adjust as needed
    print(head(processed_data[[specific_file_name]]))
}

