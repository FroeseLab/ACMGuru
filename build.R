# install.packages("roxygen2")
# install.packages("devtools")

library(roxygen2)
library(devtools)

# build documentation
# setwd("/path/to/your/package/ACMGuru")
devtools::document()

# check
devtools::check()
devtools::build()
