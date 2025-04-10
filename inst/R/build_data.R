dat_files <- dir("inst/extdata/real")
for (dat_file in dat_files) {
  dat <- read.csv(sprintf("inst/extdata/real/%s", dat_file), header = T)
  dat_nm <- gsub(".csv", "", dat_file)
  assign(dat_nm, dat)
  usethis::use_data_raw(dat_nm)
  source(sprintf("data-raw/%s.R", dat_nm))
}

