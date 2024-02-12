![Alt Text](./images/DALLE_guru.jpg)

# ACMGuru

ACMGuru is an R package designed to facilitate the interpretation of genetic determinants of disease in genomic data. It applies extensive annotation and filtering based on the ACMG (American College of Medical Genetics and Genomics) and AMP (Association for Molecular Pathology) guideline standards. 

Richards et al.'s joint consensus recommendation on the standards and guidelines for the interpretation of sequence variants is a foundational document in this field (2015, [doi: 10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)). ACMGuru is designed to implement the initial steps of our filtering protocol for the addition of ACMG-standardised labels to candidate causal variants, according to these standards.

Accepts VCF format.


## Features

- Automatic annotation and filtering of genomic data based on the ACMG/AMP guideline standards.
- Generation of standardised labels for candidate causal variants.
- Customisation options for filtering and annotation criteria.
- Comprehensive visualization tools for interpreting filtering and annotation results.

## NOTE: THIS PROJECT IS NOT COMPLETE

The following directions are for development and will be updated for users shortly.

## Installation

To install the latest development version of ACMGuru from GitHub:

```r
# install.packages("devtools")
devtools::install_github("your_github_username/ACMGuru")
```

## Usage

```
library(ACMGuru)

# Load your genomic data
genomic_data <- ...

# Apply ACMGuru annotation and filtering
result <- ACMGuru::acmg_filter(genomic_data)

# Visualize results
ACMGuru::plot_result(result)
```

## Documentation

For more detailed information on using ACMGuru, see the package vignette or the reference manual.

## Contributing

I welcome contributions of all types! Please see Contributing Guide for more details.

## Contact

Please report any issues or questions on the GitHub Issues page.

"link_to_vignette", "link_to_reference_manual", "link_to_contributing_guide", "license",
