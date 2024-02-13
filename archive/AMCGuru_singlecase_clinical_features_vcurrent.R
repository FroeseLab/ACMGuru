# AMCGuru ----
source("AMCGuru_singlecase_vcurrent.R")
# List all objects in the current environment
all_objects <- ls()

# Define the objects to keep
objects_to_keep <- c("df_report", "df_report_main_text")

# Identify objects to remove (all objects except those to keep)
objects_to_remove <- setdiff(all_objects, objects_to_keep)

# Remove the identified objects
rm(list = objects_to_remove, envir = .GlobalEnv)

# clinical ----
# The raw clinical data summary is shown in the following
setwd("../cohort_summary_curated_r/src/")
source("cohort_summary_vcurrent.R")

# merge ----
# sample.id is dropped in original script so do again:
samples <- read.csv("../data/SAMPLE_LIST", header = F)
samples$sample <- samples$V1
samples <- samples |> dplyr::select(-V1)

# Pheno ----
# Create new column "cohort_pheno"
samples$cohort_pheno <- samples$sample

# Replace any value that starts with "setpt" with "0" in the "cohort_pheno" column
samples$cohort_pheno[grep("^setpt", samples$sample)] <- "0"

# Replace any value that does not start with "setpt" with "1" in the "cohort_pheno" column
samples$cohort_pheno[!grepl("^setpt", samples$sample)] <- "1"

# clean IDs
samples <- separate(samples, sample, into = c("V1", "V2", "V3", "V4", "V5"))

samples <- samples |>
  mutate(V1 = ifelse(V1 == "raw", NA, V1))

samples <- samples |>
  unite(V1, V2, col = "sample.id", sep = "", na.rm = TRUE)

samples <- samples |> filter(cohort_pheno == 1)

# Clinical data
df <- read.csv("../../../dataset2/data/sepsis_v2.csv")
names(df)
df <- df |>dplyr::select(
  -exome_dataset_1,
  -exome_dataset_1.1,
  -exome_dataset_2_path,
  -exome_dataset_1_path,
  -sqlpkey,
  -personal.id)

df$sample.id <- gsub("-", "", df$sample.id)

df <- merge(samples, df, by = "sample.id", all.x = TRUE)

missing_samples <- subset(df, is.na(study.site))
missing_sample_ids <- missing_samples$sample.id

setwd("AMCGuru_singlecase")

# Assuming 'df' is the clinical data dataframe from 'cohort_summary_vcurrent.R'

hold1 <- df_report_main_text
df_report_main_text <- df_report_main_text %>% dplyr::select(sample.id) %>% as.data.frame()

# Merge df_report_main with clinical data based on 'sample.id'
merged_data <- merge(df_report_main_text, df, by = "sample.id", all.x = TRUE)

# write.csv(merged_data,  paste0("../../data/singlecase/AMCGuru_singlecase_df_report_main_text_clinical.csv"))

# summary ----
# merged_data is the dataframe we want to analyze

hold <- df
df <- merged_data

# Summary stats ----
# table categorical ----
# Function to calculate frequency count and percentage for a categorical variable
calculate_category_counts <- function(data, variable_name) {
  # Count the frequency of each category
  category_counts <- data %>%
    dplyr::select(all_of(variable_name)) %>%
    table() %>%
    as.data.frame()

  # Rename columns for clarity and add a variable name column
  names(category_counts) <- c("Category", "Count")
  category_counts$Variable <- variable_name

  # Calculate the percentage of each category
  total_count <- sum(category_counts$Count)
  category_counts$Percentage <- (category_counts$Count / total_count) * 100

  # Round the percentage values to two decimal places
  category_counts$Percentage <- round(category_counts$Percentage, digits = 2)

  return(category_counts)
}

# Identify character variables and exclude date variables
categorical_vars <- names(df)[sapply(df, is.character)]
date_vars <- c("hosp.adm", "hosp.dis", "picu.adm", "picu.dis", "bc.sampling", "death.date")
categorical_vars <- setdiff(categorical_vars, date_vars)

# Apply the function to all categorical variables and combine the results
all_category_counts <- lapply(categorical_vars, function(var) {
  calculate_category_counts(df, var)
}) %>% bind_rows()

all_category_counts <-
  all_category_counts %>%
  dplyr::select(Variable, Category, Count, Percentage) %>%
  filter(Category !="no")

# Order by 'Count' within each 'Variable' group
all_category_counts_ordered <- all_category_counts %>%
  group_by(Variable) %>%
  arrange(Variable, desc(Count)) %>%
  ungroup()  # Ungroup for further operations if needed


# Print the ordered combined frequency table
print(all_category_counts_ordered)

unique(all_category_counts_ordered$Variable)

# Filter out date variables
all_category_counts_ordered <- all_category_counts_ordered |>
  filter(!(Variable %in% c(  "HGVSc",
                             "HGVSp",
                             "V3",
                             "V4",
                             "V5",
                             "rownames",
                             "sample.id")))

# Save the ordered combined frequency table to a CSV file
write.csv(all_category_counts_ordered, file = "all_category_counts_singlecase.csv", row.names = FALSE)

# table continuous ----

# Convert the data from wide to long (similar to the histogram plot preparation)
df_long <- df |>
  dplyr::select(which(sapply(df, is.numeric))) |>  #dplyr::select numeric columns
  gather(key = "variable", value = "value")  # Convert from wide to long format

# Function to calculate summary statistics for each variable
calculate_continuous_summary_stats <- function(data, variable_name) {
  variable_data <- data %>% filter(variable == variable_name)

  summary_stats <- variable_data %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Min = min(value, na.rm = TRUE),
      Max = max(value, na.rm = TRUE),
      N = n()
    )

  # Round the statistics to two decimal places
  summary_stats$Mean <- round(summary_stats$Mean, digits = 2)
  summary_stats$Median <- round(summary_stats$Median, digits = 2)
  summary_stats$SD <- round(summary_stats$SD, digits = 2)
  summary_stats$Min <- round(summary_stats$Min, digits = 2)
  summary_stats$Max <- round(summary_stats$Max, digits = 2)


  # Add a variable name column
  summary_stats$Variable <- variable_name

  return(summary_stats)
}

# Apply the function to each unique variable in the df_long dataframe
all_continuous_summary_stats <- lapply(unique(df_long$variable), function(var) {
  calculate_continuous_summary_stats(df_long, var)
}) %>% bind_rows()

# Rearrange columns for better readability
all_continuous_summary_stats <- all_continuous_summary_stats %>%
  dplyr::select(Variable, Mean, Median, SD, Min, Max, N)

# Print the combined summary statistics table
print(all_continuous_summary_stats)

# Save the combined summary statistics table to a CSV file
write.csv(all_continuous_summary_stats, file = "all_continuous_stats_singlecase.csv", row.names = FALSE)

