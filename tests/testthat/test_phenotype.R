# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- test phenotype class.


# set.seed(2110) # for reproductible RNG


##### Initialisation functions ####
if (interactive()) {
  devtools::load_all()
  source("tests/testthat/src/functionsForTests.R")
} else source("src/functionsForTests.R")




#### TESTS  traits ####
test_that("quant trait initialisation", {
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec, nMarker = 300)

  expect_error({
    myTrait <- trait$new(
      name = "myTrait",
      qtn = sample(SNPs$SNPcoord$SNPid, 100),
      qtnEff = rnorm(100, sd = 0.5)
    )
  }, NA)
  expect_error({
    qtn <- sample(SNPs$SNPcoord$SNPid, 100)
    qtnEff <- rnorm(100, sd = 0.5)
    names(qtnEff) <- sample(qtn, length(qtn))
    myTrait2 <- trait$new(
      name = "myTrait",
      qtn = qtn,
      qtnEff = qtnEff
    )
  }, NA)
  expect_equal(names(myTrait$qtnEff), myTrait$qtn)
  expect_equal(names(myTrait2$qtnEff), myTrait2$qtn)

})


test_that("quant trait gv calculation", {
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec, nMarker = 300)
  nInds <- 100
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })
  indList <- create_inds(haploList)
  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)
  myTrait <- trait$new(
    name = "myTrait",
    qtn = sample(SNPs$SNPcoord$SNPid, 100),
    qtnEff = rnorm(100, sd = 0.5)
  )

  expect_error({gv <- myTrait$gv(myPop)}, NA)
  expect_is(gv, "matrix")
  expect_equal(nrow(gv), nInds)
  expect_equal(row.names(gv), names(myPop$inds))

})





#### TESTS  phenotyper ####
test_that("phenotyper initialisation", {
  mySpec <- create_spec()
  SNPs <- create_SNP(mySpec, nMarker = 300)
  nInds <- 100
  haploList <- lapply(seq(nInds), function(x){
    create_haplo(SNPs)
  })

  indList <- create_inds(haploList)

  myPop <- population$new(name = "My Population 1",
                          inds = indList,
                          verbose = FALSE)

  myTrait1 <- trait$new(
    name = "myTrait1",
    qtn = sample(SNPs$SNPcoord$SNPid, 200),
    qtnEff = rnorm(200, sd = 0.5)
  )
  myTrait2 <- trait$new(
    name = "myTrait2",
    qtn = sample(SNPs$SNPcoord$SNPid, 100),
    qtnEff = rnorm(100, sd = 0.2)
  )
  myTrait3 <- trait$new(
    name = "myTrait3",
    qtn = sample(SNPs$SNPcoord$SNPid, 50),
    qtnEff = rnorm(50, sd = 0.05)
  )
  myTraits <- list(myTrait1, myTrait2, myTrait3)

  expect_error({
    phenoLab1 <- phenotyper$new(name = "My phenoLab",
                               traits = myTraits,
                               plotCost = 1,
                               mu = c(0, 10, 50),
                               ve = NULL,
                               he = c(0.4, 0.5, 0.6),
                               pop = myPop)
  }, NA)

  expect_error({
    phenoLab2 <- phenotyper$new(name = "My phenoLab",
                               traits = myTraits,
                               plotCost = 1,
                               mu = c(0, 10, 50),
                               ve = c(20, 5, 1),
                               he = NULL,
                               pop = NULL)
  }, NA)

  expect_error({
    phenoLab1$he(myPop)
  }, NA)
  expect_error({
    phenoLab2$he(myPop)
  }, NA)

  expect_error({
    pheno <- phenoLab1$trial(myPop, rep = 2, offset = 4)
  }, NA)
  expect_error({
    phenoLab2$trial(myPop, rep = 2, offset = c(-3, -4, 5))
  }, NA)
  expect_error({
    phenoLab1$trial(myPop, rep = round(runif(myPop$nInd, 1, 3)), offset = 4)
  })


  expect_equal(names(pheno$data),
               c("ind", "myTrait1", "myTrait2", "myTrait3", "rep", "phenotyper"))
  expect_true(!any(is.na(pheno$data$myTrait1)))


})