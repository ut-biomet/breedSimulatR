# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- test individual crossing.


# set.seed(8065) # for reproductible RNG


##### Initialisation functions ####
if (basename(getwd()) == "breedSimulatR") {
  devtools::load_all()
  source("tests/testthat/src/functionsForTests.R")
} else source("src/functionsForTests.R")


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
  expect_error(makeSingleCross(myInd1, myInd2, "newInd"), NA)
  expect_error(makeSingleCross(myInd1, myInd2, "newInd", 2), NA)

  expect_is(makeSingleCross(myInd1, myInd2, "newInd"), "list")
  expect_equal(length(makeSingleCross(myInd1, myInd2, "newInd")), 1)
  expect_equal(length(makeSingleCross(myInd1, myInd2, "newInd", n = 3)), 3)
  expect_is(makeSingleCross(myInd1, myInd2, "newInd")[[1]], "individual")

  expect_error(makeSingleCross(myInd1, myInd2))

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
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)

  nCross <- 5
  crossToDo <- data.frame(ind1 = sample(names(myPop$inds), nCross, replace = T),
                          ind2 = sample(names(myPop$inds), nCross, replace = T),
                          n = rep_len(c(1,2,3), length.out = nCross),
                          names = paste("newInd", 1:nCross))


  #### Tests:
  expect_error({newInds <- makeCrosses(crossToDo, myPop)}, NA)
  expect_is(newInds, "list")
  expect_equal(length(newInds), sum(crossToDo$n))
  lapply(newInds, function(x){
    expect_is(x, "individual")
  })


  # check offspring names are unique
  offNames <- as.character(vapply(newInds, function(x){x$name}, "character"))
  expect_is(offNames, "character")
  expect_true(length(offNames) == length(unique(offNames)))

})

test_that("makeCrosses particular cases", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  nInds <- 10
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)


  #### Tests:
  # 1 offspring per cross
  nCross <- 3
  crossToDo <- data.frame(ind1 = sample(names(myPop$inds), nCross),
                          ind2 = sample(names(myPop$inds), nCross),
                          n = 1,
                          names = NA)

  expect_error({listOff <- makeCrosses(crossToDo, myPop)}, NA)
  expect_is(listOff, "list")
  lapply(listOff, function(x){
    expect_is(x, "individual")
  })

  # all offsprings with differents names when they all have the same parents
  nCross <- 3
  crossToDo <- data.frame(ind1 = rep(sample(names(myPop$inds), 1), nCross),
                          ind2 = rep(sample(names(myPop$inds), 1), nCross),
                          n = c(1,1,5),
                          names = NA)
  listOff <- makeCrosses(crossToDo, myPop)
  offNames <- as.character(vapply(listOff, function(x){x$name}, "character"))
  expect_true(length(offNames) == length(unique(offNames)))

})

test_that("makeCrosses errors, messages", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  nInds <- 10
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)
  nCross <- 5
  crossToDo <- data.frame(ind1 = sample(names(myPop$inds), nCross, replace = T),
                          ind2 = sample(names(myPop$inds), nCross, replace = T),
                          n = rep_len(c(1,2,3), length.out = nCross),
                          names = paste("newInd", 1:nCross),
                          stringsAsFactors = FALSE)
  crossToDoInit <- crossToDo


  #### Tests
  # test NA in crossToDo$ind1
  crossToDo[sample(nCross, 1), "ind1"] <- NA
  expect_error(makeCrosses(crossToDo, myPop),
               'Columns "ind1" and "ind2" should not contain any "NA"')
  crossToDo <- crossToDoInit
  # test NA in crossToDo$ind2
  crossToDo[sample(nCross, 1), "ind2"] <- NA
  expect_error(makeCrosses(crossToDo, myPop),
               'Columns "ind1" and "ind2" should not contain any "NA"')
  # test NA in crossToDo$names
  crossToDo[sample(nCross, 1), "names"] <- NA
  expect_error(makeCrosses(crossToDo, myPop),
               'Columns "ind1" and "ind2" should not contain any "NA"')
  crossToDo <- crossToDoInit



  # test parents ind1 not in population
  crossToDo[sample(nCross, 1), "ind1"] <- "toto"
  expect_error(makeCrosses(crossToDo, myPop),
               'Parents not found in the population')
  crossToDo <- crossToDoInit
  # test parents ind2 not in population
  crossToDo[sample(nCross, 1), "ind2"] <- "tata"
  expect_error(makeCrosses(crossToDo, myPop),
               'Parents not found in the population')
  crossToDo <- crossToDoInit


  # offsprings 'names in population
  crossToDo$names <- crossToDo$ind1
  expect_message(makeCrosses(crossToDo, myPop),
               'Offspring names already exist in the population')


})
