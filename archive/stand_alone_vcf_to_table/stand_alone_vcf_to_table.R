# Set paramenters
# file_suffix <- "output4"

# require(dplyr)
# library(tidyr)

# /////////////////////////////////////////////////////////////////////////////
# Set filter parameters ----
# /////////////////////////////////////////////////////////////////////////////

# VCF file has 940 samples. 
CRITICAL_n_samples <- 940
gnomad_total <- 140000 
gnomad_rare_thresh <- 1 # keep all
gnomad_count_thresh <- gnomad_total*gnomad_rare_thresh
cohort_max_carriers <- 940 # keep all
cohort_freq_max <- cohort_max_carriers/CRITICAL_n_samples

# /////////////////////////////////////////////////////////////////////////////
# Run VCF to table scripts ----
# /////////////////////////////////////////////////////////////////////////////

# for (f in 1:length(file_list)) {
	source("./stand_alone_vcf_to_table/gather.R")
	source("./stand_alone_vcf_to_table/genotype_clean.R")
	source("./stand_alone_vcf_to_table/progress_bar.R")
	rm(df_main) # clean if re-running
# }

# names(df)
