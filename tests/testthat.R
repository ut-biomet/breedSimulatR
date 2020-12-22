library(testthat)
library(breedSimulatR)

options(testthat.progress.max_fails = 200)
test_check("breedSimulatR", stop_on_warning = FALSE)
