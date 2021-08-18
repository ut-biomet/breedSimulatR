# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- test individual class.


# set.seed(2110) # for reproductible RNG


##### Initialisation functions ####
if (basename(getwd()) == "breedSimulatR") {
  devtools::load_all()
  source("tests/testthat/src/functionsForTests.R")
} else source("src/functionsForTests.R")



#### TESTS ####
test_that("individual initialisation", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  haplo <- create_haplo(SNPs)

  # create individual
  expect_error({myInd <-  individual$new(name = "Ind 1",
                              specie = mySpec,
                              parent1 = "OkaaSan",
                              parent2 = "OtouSan",
                              haplo = haplo,
                              verbose = F)},
               NA)

  expect_identical(myInd$name, "Ind 1")
  expect_identical(myInd$specie, mySpec)
  expect_identical(myInd$parent1, "OkaaSan")
  expect_identical(myInd$parent2, "OtouSan")
  expect_identical(myInd$haplo, haplo)
})


test_that("individual gametes generation", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  haplo <- create_haplo(SNPs)

  # create individual
  myInd <- individual$new(name = "Ind 1",
                          specie = mySpec,
                          parent1 = "OkaaSan",
                          parent2 = "OtouSan",
                          haplo = haplo,
                          verbose = F)

  #### Tests
  expect_error(myInd$generateGametes(), NA)

  # check output class
  expect_is(myInd$generateGametes(), "list")
  expect_equal(length(myInd$generateGametes()), 1)
  expect_equal(length(myInd$generateGametes(3)), 3)

  # check each gamete
  expect_is(myInd$generateGametes()[[1]], "integer")
  expect_equal(length(myInd$generateGametes()[[1]]), SNPs$nSNP())

  # check for markers names
  expect_equal(names(myInd$generateGametes()[[1]]),
               names(myInd$haplo$allelDose))

  # check gamete genotype corresponds to parent's genotype:
  gam <- myInd$generateGametes()[[1]]
  c1 <- gam == do.call(cbind, myInd$haplo$values)[1,]
  c2 <- gam == do.call(cbind, myInd$haplo$values)[2,]
  expect_true(all(c1 | c2))
})
