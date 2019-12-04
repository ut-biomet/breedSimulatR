# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Test File -- test haplotype class.


# set.seed(8241) # for reproductible RNG


test_that("haplotype initialisation", {
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

  # create haplotype object
  expect_error({haplo <- haplotype$new(SNPinfo = SNPs,
                                       haplo = rawHaplo)},
               NA)

  expect_is(haplo$SNPinfo, "SNPinfo")
  expect_is(haplo$values, "list")
  expect_equal(length(haplo$values), 3)
  expect_equal(names(haplo$values), c("Chr1", "Chr2", "Chr3"))
  expect_is(unlist(haplo$values), "integer")
  expect_equal(haplo$allelDose, colSums(rawHaplo)[names(haplo$allelDose)])
})

test_that("haplotype initialisation errors", {
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



  expect_error(haplotype$new(SNPinfo = SNPs,
                             haplo = rawHaplo),
               "haplo must be a named matrix")

  rawHaplo <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = T), nrow = 2)
  colnames(rawHaplo) <- rep("toto", ncol(rawHaplo))
  expect_error(haplotype$new(SNPinfo = SNPs,
                             haplo = rawHaplo),
               "colnames\\(haplo\\) must be the names of the markers")
})
