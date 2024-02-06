

# varsome ----
# LE = less than equal to, GE = greater than equal to
# varsome <- read.csv(file = "../../data/singlecase/varsome_calibrated_insilico_thresholds.tsv", sep="\t")


# source("../R/acmguru_vcurrent.R")

if (!dir.exists("./output/")) {
  dir.create("./output/", recursive = TRUE)
}

# Specify the path for input data and samples file
input_path <- "./data/"
# input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"
samples_file_path <- "./data/samples.tsv"
af_threshold <- 0.1  # Allele frequency threshold
gnomad_max <- 1e-6

# In your test.R script
library(devtools)
load_all("../R/") # Adjust the path to the root of your package


# check metadataclass_file_path

# Merge all dataframes in the list into a single dataframe
all_data <- bind_rows(processed_data_list)
#
# processed_data_list[[1]] |> nrow()
# processed_data_list[[2]] |> nrow()
# processed_data_list[[1]] |> ncol()
# processed_data_list[[2]] |> ncol()
#
# processed_data_list[[1]]$BayesDel_addAF_pred |> class()
# processed_data_list[[2]]$BayesDel_addAF_pred |> class()
#
# processed_data_list[[1]]$BayesDel_addAF_pred |> head()
# processed_data_list[[2]]$BayesDel_addAF_pred |> head()
#
# df1 <- processed_data_list[[1]]
# df2 <- processed_data_list[[2]]
#
# processed_data_list[[1]]$MutationTaster_score |> class()
# processed_data_list[[2]]$MutationTaster_score |> class()
#
# processed_data_list[[1]]$MutationTaster_score |> head()
# processed_data_list[[2]]$MutationTaster_score |> head()


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
