# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- test individual crossing.


# set.seed(8065) # for reproductible RNG


##### Initialisation functions ####
source("src/functionsForTests.R")



#### TESTS ####
test_that("makeSingleCross", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  haplo1 <- create_haplo(SNPs)
  haplo2 <- create_haplo(SNPs)
  myInds <- create_inds(list(haplo1, haplo2))
  myInd1 <- myInds[[1]]
  myInd2 <- myInds[[2]]



  #### Tests:
  expect_error(makeSingleCross(myInd1, myInd2), NA)
  expect_is(makeSingleCross(myInd1, myInd2), "list")
  expect_equal(length(makeSingleCross(myInd1, myInd2)), 1)
  expect_equal(length(makeSingleCross(myInd1, myInd2, n = 3)), 3)
  expect_is(makeSingleCross(myInd1, myInd2)[[1]], "individual")

  newInd <- makeSingleCross(myInd1, myInd2, names = "New", verbose = F)[[1]]
  expect_identical(newInd$name, "New")
  expect_identical(newInd$specie, mySpec)
  expect_identical(newInd$parent1, "Ind 1")
  expect_identical(newInd$parent2, "Ind 2")

  # check genotype:
  mustEq2 <- myInd1$haplo$allelDose == 2 & myInd2$haplo$allelDose == 2
  mustEq0 <- myInd1$haplo$allelDose == 0 & myInd2$haplo$allelDose == 0
  mustEq1.a <- myInd1$haplo$allelDose == 2 & myInd2$haplo$allelDose == 0
  mustEq1.b <- myInd1$haplo$allelDose == 0 & myInd2$haplo$allelDose == 2
  mustEq12.a <- myInd1$haplo$allelDose == 2 & myInd2$haplo$allelDose == 1
  mustEq12.b <- myInd1$haplo$allelDose == 1 & myInd2$haplo$allelDose == 2
  mustEq01.a <- myInd1$haplo$allelDose == 0 & myInd2$haplo$allelDose == 1
  mustEq01.b <- myInd1$haplo$allelDose == 1 & myInd2$haplo$allelDose == 0
  mustEq012 <- myInd1$haplo$allelDose == 1 & myInd2$haplo$allelDose == 1

  expect_true(all(newInd$haplo$allelDose[mustEq2] == 2))
  expect_true(all(newInd$haplo$allelDose[mustEq0] == 0))
  expect_true(all(newInd$haplo$allelDose[mustEq1.a] == 1))
  expect_true(all(newInd$haplo$allelDose[mustEq1.b] == 1))
  expect_true(all(newInd$haplo$allelDose[mustEq12.a] %in% c(1,2)))
  expect_true(all(newInd$haplo$allelDose[mustEq12.b] %in% c(1,2)))
  expect_true(all(newInd$haplo$allelDose[mustEq01.a] %in% c(1,0)))
  expect_true(all(newInd$haplo$allelDose[mustEq01.b] %in% c(1,0)))
  expect_true(all(newInd$haplo$allelDose[mustEq012] %in% c(0,1,2)))

})


test_that("makeCrosses", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  nInds <- 10
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)


  # nOff <- 2
  # crossToDo <- data.frame(ind1 = sample(IndNames, nOff, replace = T),
  #                         ind2 = sample(IndNames, nOff, replace = T),
  #                         n = 1,
  #                         names = paste("newInd", 1:nOff))

  # newPop <- makeCrosses(crossToDo, pop)
  #
  # testOffspring <- function(newInd){
  #   parent1 <- pop[[newInd$parent1]]
  #   parent2 <- pop[[newInd$parent2]]
  #
  #   expect_identical(newInd$specie, mySpec)
  #   expect_identical(newInd$parent1, parent1$name)
  #   expect_identical(newInd$parent2, parent2$name)
  #
  #   # check genotype:
  #   mustEq2 <- parent1$haplo$allelDose == 2 & parent2$haplo$allelDose == 2
  #   mustEq0 <- parent1$haplo$allelDose == 0 & parent2$haplo$allelDose == 0
  #   mustEq1.a <- parent1$haplo$allelDose == 2 & parent2$haplo$allelDose == 0
  #   mustEq1.b <- parent1$haplo$allelDose == 0 & parent2$haplo$allelDose == 2
  #   mustEq12.a <- parent1$haplo$allelDose == 2 & parent2$haplo$allelDose == 1
  #   mustEq12.b <- parent1$haplo$allelDose == 1 & parent2$haplo$allelDose == 2
  #   mustEq01.a <- parent1$haplo$allelDose == 0 & parent2$haplo$allelDose == 1
  #   mustEq01.b <- parent1$haplo$allelDose == 1 & parent2$haplo$allelDose == 0
  #   mustEq012 <- parent1$haplo$allelDose == 1 & parent2$haplo$allelDose == 1
  #
  #   expect_true(all(newInd$haplo$allelDose[mustEq2] == 2))
  #   expect_true(all(newInd$haplo$allelDose[mustEq0] == 0))
  #   expect_true(all(newInd$haplo$allelDose[mustEq1.a] == 1))
  #   expect_true(all(newInd$haplo$allelDose[mustEq1.b] == 1))
  #   expect_true(all(newInd$haplo$allelDose[mustEq12.a] %in% c(1,2)))
  #   expect_true(all(newInd$haplo$allelDose[mustEq12.b] %in% c(1,2)))
  #   expect_true(all(newInd$haplo$allelDose[mustEq01.a] %in% c(1,0)))
  #   expect_true(all(newInd$haplo$allelDose[mustEq01.b] %in% c(1,0)))
  #   expect_true(all(newInd$haplo$allelDose[mustEq012] %in% c(0,1,2)))
  # }
  #
  #
  # #### Tests:
  # lapply(newPop, testOffspring)

})