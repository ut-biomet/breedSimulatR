# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- test individual crossing.


# set.seed(8065) # for reproductible RNG


test_that("makeSingleCross", {

  #### Initialisation:
  # create specie
  nChr <- round(runif(1, 1, 10))
  lchr <- round(rnorm(nChr, 450, 50))
  ploidy <- 2
  recombRate <- 3/sum(lchr) # expect 3 recombination event

  mySpec <- specie$new(nChr = nChr,
                       lchr = lchr,
                       ploidy = ploidy,
                       recombRate = recombRate,
                       verbose = F)

  # simulate SNP
  nMarker <- round(lchr/10)

  # generate arbitrary marker position
  pos <- c()
  for (i in seq(nChr)) {
    pos <- c(pos, sample(lchr[i], nMarker[i]))
  }

  # generate arbitrary SNPid
  SNPid <- sprintf(fmt = paste0("SNP%0", ceiling(log10(sum(nMarker))),"i"),
                   sample(5000, sum(nMarker)))
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarker),
                         pos = pos,
                         SNPid = SNPid)
  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  # simulate haplotypes
  rawHaplo1 <- matrix(sample(c(0, 1), sum(nMarker) * ploidy, replace = T),
                     nrow = ploidy)
  colnames(rawHaplo1) <- SNPid
  haplo1 <- haplotype$new(SNPinfo = SNPs,
                         haplo = rawHaplo1)
  rawHaplo2 <- matrix(sample(c(0, 1), sum(nMarker) * ploidy, replace = T),
                      nrow = ploidy)
  colnames(rawHaplo2) <- SNPid
  haplo2 <- haplotype$new(SNPinfo = SNPs,
                          haplo = rawHaplo2)

  # create individuals
  myInd1 <- individual$new(name = "Ind 1",
                          specie = mySpec,
                          parent1 = "OkaaSan1",
                          parent2 = "OtouSan1",
                          haplo = haplo1,
                          verbose = F)
  myInd2 <- individual$new(name = "Ind 2",
                           specie = mySpec,
                           parent1 = "OkaaSan2",
                           parent2 = "OtouSan2",
                           haplo = haplo2,
                           verbose = F)


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

  #### Initialisation:
  # create specie
  nChr <- round(runif(1, 1, 10))
  lchr <- round(rnorm(nChr, 450, 50))
  ploidy <- 2
  recombRate <- 3/sum(lchr) # expect 3 recombination event

  mySpec <- specie$new(nChr = nChr,
                       lchr = lchr,
                       ploidy = ploidy,
                       recombRate = recombRate,
                       verbose = F)

  # simulate SNP
  nMarker <- round(lchr/10)

  # generate arbitrary marker position
  pos <- c()
  for (i in seq(nChr)) {
    pos <- c(pos, sample(lchr[i], nMarker[i]))
  }

  # generate arbitrary SNPid
  SNPid <- sprintf(fmt = paste0("SNP%0", ceiling(log10(sum(nMarker))),"i"),
                   sample(5000, sum(nMarker)))
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarker),
                         pos = pos,
                         SNPid = SNPid)
  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  # simulate haplotypes
  nInds <- 10
  rawHaploList <- lapply(c(1:nInds), function(x){
    rawHaplo <- matrix(sample(c(0, 1), sum(nMarker) * ploidy, replace = T),
                       nrow = ploidy)
    colnames(rawHaplo) <- SNPid
    rawHaplo
  })

  haploList <- lapply(rawHaploList, function(rawHaplo){
    haplo <- haplotype$new(SNPinfo = SNPs,
                           haplo = rawHaplo)
  })

  pop <- mapply(function(haplo, id){
    myInd <- individual$new(name = paste("Ind", id),
                            specie = mySpec,
                            parent1 = paste("OkaaSan", id),
                            parent2 = paste("OtouSan", id),
                            haplo = haplo,
                            verbose = F)
  },
  haploList, c(1:nInds))
  names(pop) <- paste("Ind", c(1:nInds))

  nOff <- 2
  crossToDo <- data.frame(ind1 = sample(names(pop), nOff, replace = T),
                          ind2 = sample(names(pop), nOff, replace = T),
                          n = 1,
                          names = paste("newInd", 1:nOff))

  newPop <- makeCrosses(crossToDo, pop)

  testOffspring <- function(newInd){
    parent1 <- pop[[newInd$parent1]]
    parent2 <- pop[[newInd$parent2]]

    expect_identical(newInd$specie, mySpec)
    expect_identical(newInd$parent1, parent1$name)
    expect_identical(newInd$parent2, parent2$name)

    # check genotype:
    mustEq2 <- parent1$haplo$allelDose == 2 & parent2$haplo$allelDose == 2
    mustEq0 <- parent1$haplo$allelDose == 0 & parent2$haplo$allelDose == 0
    mustEq1.a <- parent1$haplo$allelDose == 2 & parent2$haplo$allelDose == 0
    mustEq1.b <- parent1$haplo$allelDose == 0 & parent2$haplo$allelDose == 2
    mustEq12.a <- parent1$haplo$allelDose == 2 & parent2$haplo$allelDose == 1
    mustEq12.b <- parent1$haplo$allelDose == 1 & parent2$haplo$allelDose == 2
    mustEq01.a <- parent1$haplo$allelDose == 0 & parent2$haplo$allelDose == 1
    mustEq01.b <- parent1$haplo$allelDose == 1 & parent2$haplo$allelDose == 0
    mustEq012 <- parent1$haplo$allelDose == 1 & parent2$haplo$allelDose == 1

    expect_true(all(newInd$haplo$allelDose[mustEq2] == 2))
    expect_true(all(newInd$haplo$allelDose[mustEq0] == 0))
    expect_true(all(newInd$haplo$allelDose[mustEq1.a] == 1))
    expect_true(all(newInd$haplo$allelDose[mustEq1.b] == 1))
    expect_true(all(newInd$haplo$allelDose[mustEq12.a] %in% c(1,2)))
    expect_true(all(newInd$haplo$allelDose[mustEq12.b] %in% c(1,2)))
    expect_true(all(newInd$haplo$allelDose[mustEq01.a] %in% c(1,0)))
    expect_true(all(newInd$haplo$allelDose[mustEq01.b] %in% c(1,0)))
    expect_true(all(newInd$haplo$allelDose[mustEq012] %in% c(0,1,2)))
  }


  #### Tests:
  lapply(newPop, testOffspring)

})