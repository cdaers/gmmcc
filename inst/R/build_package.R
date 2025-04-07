# build the help files
roxygen2::roxygenise(clean = TRUE)

# build data: run build_data.R
source("inst/R/build_data.R")

# build package: Build -> Build Binary Package
