# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- Selection functions.


# set.seed(7352) # for reproductible RNG


##### Initialisation functions ####
if (interactive()) {
  devtools::load_all()
  source("tests/testthat/src/functionsForTests.R")
} else source("src/functionsForTests.R")



#### TESTS selectBV ####
test_that("selectBV", {
  #### Initialisation
  mySpec <- create_spec(lchr = 5000)
  SNPs <- create_SNP(mySpec)
  nInds <- 100
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)

  snpEffects <- rnorm(SNPs$nSNP(),mean = 0, sd = 42/sqrt(SNPs$nSNP()))
  names(snpEffects) <- colnames(myPop$genoMat)

  nSel <- 10

  #### Tests:
  expect_error({selectedInds <- selectBV(myPop, nSel, snpEffects)},
               NA)
  expect_is(selectedInds, "character")
  expect_equal(length(selectedInds), nSel)
  expect_equal(unique(selectedInds), selectedInds)

  bv <- myPop$genoMat %*% snpEffects
  bv <- as.numeric(bv)
  names(bv) <- rownames(myPop$genoMat)

  expect_equivalent(bv[selectedInds[1]], max(bv))

  bv <- sort(bv, decreasing = T)
  expect_equal(selectedInds, names(bv)[1:nSel])

})


#### TESTS selectWBV ####
test_that("selectWBV", {
  #### Initialisation
  mySpec <- create_spec(lchr = 5000)
  SNPs <- create_SNP(mySpec)
  nInds <- 100
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)

  snpEffects <- rnorm(SNPs$nSNP(),mean = 0, sd = 42/sqrt(SNPs$nSNP()))
  names(snpEffects) <- colnames(myPop$genoMat)

  nSel <- 10

  #### Tests:
  expect_error({selectedInds <- selectWBV(myPop, nSel, snpEffects)},
               NA)
  expect_is(selectedInds, "character")
  expect_equal(length(selectedInds), nSel)
  expect_equal(unique(selectedInds), selectedInds)

  warning("calculation's results not checked !")

})