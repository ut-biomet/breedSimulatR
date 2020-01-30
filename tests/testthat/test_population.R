# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- test population class.


# set.seed(5754) # for reproductible RNG


##### Initialisation functions ####
source("src/functionsForTests.R")



#### TESTS ####
test_that("population initialisation", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  nInds <- 5
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)


  #### Tests:
  expect_error({myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)},
               NA)
  expect_equal(myPop$nInd, nInds)
  expect_is(myPop$inds, "list")
  expect_equal(length(myPop$inds), myPop$nInd)
  expect_is(myPop$genoMat, "matrix")
  expect_equal(dim(myPop$genoMat), c(myPop$nInd, SNPs$nSNP()))

  # check individuals' id in the object pop (names(myPop$inds)) do not depend of
  # the id of the individuals in the input list (names(indList)) :
  names(indList) <- paste("toto", 1:length(indList))
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)
  expect_equal(names(myPop$inds),
               as.character(sapply(indList, function(x){x$name})))


})



test_that("population initialisation particular cases", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  nInds <- 5
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)

  # check initialisation with one individual
  expect_error({myPop <- population$new(name = "My Population 1",
                                        inds = indList[[1]],
                                        verbose = FALSE)},
               NA)


})

test_that("population add individuals", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)

  haplo1 <- create_haplo(SNPs)
  haplo2 <- create_haplo(SNPs)
  ind1 <-  create_inds(haplo1)
  ind2 <-  create_inds(haplo2)
  ind2$name <- "Ind 2"


  #### Tests:
  expect_error({myPop <- population$new(name = "My Population 1",
                                        inds = list(ind1, ind2),
                                        verbose = FALSE)},
               NA)

  ## Check individuals with same name
  ind2$name <- "Ind 1"
  expect_error({myPop <- population$new(name = "My Population 1",
                                        inds = list(ind1, ind2),
                                        verbose = FALSE)},
                "Individual with the same name already exists in the population:")



  ## Check individuals of differents species
  mySpec1 <- create_spec(name = "Spec 1")
  mySpec2 <- create_spec(name = "Spec 2")
  SNPs1 <- create_SNP(mySpec1)
  SNPs2 <- create_SNP(mySpec2)
  haplo1 <- create_haplo(SNPs1)
  haplo2 <- create_haplo(SNPs2)
  ind1 <-  create_inds(haplo1)
  ind2 <-  create_inds(haplo2)
  ind2$name <- "Ind 2"

  expect_error({myPop <- population$new(name = "My Population 1",
                                        inds = list(ind1, ind2),
                                        verbose = FALSE)},
               "different species")
})



#### TESTS ####
test_that("population remove individuals", {
  #### Initialisation
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec)
  nInds <- 5
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <-  create_inds(haploList)

  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)

  #### Tests:
  expect_error(myPop$remInds("Ind 3"), NA)
  expect_warning(myPop$remInds("Ind 3"),
                 "Some individuals to remove are not in the population:")
  expect_error(myPop$remInds(names(myPop$inds)), NA)
  expect_warning(myPop$remInds("Ind 1"),
                 "Some individuals to remove are not in the population:")

})


test_that("population creation", {
  #### Initialisation
  mySpec <- create_spec(nChr = 10, lchr = 10^6)

  if (interactive()) {
    snpCoord <- read.csv(file = "tests/testthat/src/snpCoord.csv", header = T)
  } else snpCoord <- read.csv(file = "src/snpCoord.csv", header = T)
  SNPs <- SNPinfo$new(SNPcoord = snpCoord, specie = mySpec)

  if (interactive()) {
    geno <- read.csv(file = "tests/testthat/src/genotype.csv",
                     header = T,
                     row.names = 1)
  } else {
    geno <- read.csv(file = "src/genotype.csv",
                     header = T,
                     row.names = 1)
  }
  expect_true(all(dim(geno) > 2))
  expect_equal(ncol(geno), SNPs$nSNP())


  # TEST:
  expect_error({myPop <- createPop(geno = geno,
                                   SNPinfo = SNPs,
                                   indNames = NULL,
                                   popName = "My pop",
                                   verbose = FALSE)},
               NA)
  expect_equal(myPop$nInd, nrow(geno))
  expect_equal(names(myPop$inds), row.names(geno))
  expect_equal(myPop$genoMat, as.matrix(geno[, sort(colnames(geno))]))

})