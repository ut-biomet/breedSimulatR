# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Test File -- test individual class.


# set.seed(2110) # for reproductible RNG

test_that("individual initialisation", {
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       ploidy = 2,
                       recombRate = 0.004,
                       verbose = F)

  # simulate SNP
  SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
                                 rep("Chr2", 4),
                                 rep("Chr3", 5)),
                         pos = c(c(1,sample(100, 2)),
                                 c(3,sample(150, 3)),
                                 c(2,sample(200, 4))),
                         SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
                                         1:(3 + 4 + 5)))
  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  # simulate haplotype
  rawHaplo <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = T), nrow = 2)
  colnames(rawHaplo) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
                                1:(3 + 4 + 5))
  haplo <- haplotype$new(SNPinfo = SNPs,
                         haplo = rawHaplo)

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

  # simulate haplotype
  rawHaplo <- matrix(sample(c(0, 1), sum(nMarker) * ploidy, replace = T),
                     nrow = ploidy)
  colnames(rawHaplo) <- SNPid
  haplo <- haplotype$new(SNPinfo = SNPs,
                         haplo = rawHaplo)

  # create individual
  myInd <- individual$new(name = "Ind 1",
                          specie = mySpec,
                          parent1 = "OkaaSan",
                          parent2 = "OtouSan",
                          haplo = haplo,
                          verbose = F)

  expect_is(myInd$generateGametes(), "list")
  expect_equal(length(myInd$generateGametes()), 1)
  expect_equal(length(myInd$generateGametes(3)), 3)
  expect_equal(length(unlist(myInd$generateGametes())), SNPs$nSNP())

  # check for markers names of the gamete
  for (chr in mySpec$chrNames) {
    expect_equal(names(myInd$generateGametes()[[1]][[chr]]),
                 colnames(myInd$haplo$values[[chr]]))
  }

  geno <- unlist(myInd$generateGametes())
  names(geno) <- sub(".*(?=\\.).", "", names(geno), perl = TRUE)

  # test the markers names are similar (nessesary for the next test)
  expect_true(all(colnames(do.call(cbind, myInd$haplo$values)) == names(geno)))
  c1 <- geno == do.call(cbind, myInd$haplo$values)[1,]
  c2 <- geno == do.call(cbind, myInd$haplo$values)[2,]

  # check gamete genotype is in the individual genotye
  expect_true(all(c1 | c2))

})