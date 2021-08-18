library(testthat)
library(breedSimulatR)

options(testthat.progress.max_fails = 1000)
test_check("breedSimulatR", stop_on_warning = FALSE)


# options(testthat.progress.max_fails = 1)
# i <- 1
# seed <- 1234
# while(T){
#   seed <- round(runif(1,1,10000))
#   print(paste(i, "seed:", seed))
#   set.seed(seed)
#   expect_error({
#     x <- testthat::test_local("/home/julien/Documents/Work/1_Projects/BreedingSimulation/breedSimulatR",
#                               filter = NULL,
#                               reporter = c("stop"),
#                               stop_on_failure = TRUE,
#                               stop_on_warning = TRUE)
#   }, NA)
#   i <- i+1
# }
