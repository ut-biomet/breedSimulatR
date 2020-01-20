# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test file for the SNPinfo class


# set.seed(6705) # for reproductible RNG


test_that("SNPinfo initialization", {
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       verbose = F)

  # simulate SNP
  SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
                                 rep("Chr2", 4),
                                 rep("Chr3", 5)),
                         pos = c(c(1,sample(100, 2)),
                                 c(3,sample(150, 3)),
                                 c(2,sample(200, 4))),
                         SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
                                         1:(3 + 4 + 5)),
                         stringsAsFactors = FALSE)

  # create SNPinfo object
  expect_error({SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)},
               NA)

  expect_identical(SNPs$specie, mySpec)
  expect_is(SNPs$SNPcoord, "data.frame")
  expect_is(SNPs$SNPcoordList, "list")
  expect_equal(length(SNPs$SNPcoordList), 3)
  expect_identical(names(SNPs$SNPcoordList), mySpec$chrNames)
  expect_equal(SNPs$nSNP(), 3 + 4 + 5)

  # check SNPs$SNPcoord is sorted increasingly for all chromosomes (needed for
  # the function "findInterval" in individual's generateGametes method)
  for (chr in mySpec$chrNames) {
    expect_true(!is.unsorted(SNPs$SNPcoord[SNPs$SNPcoord$chr == chr, "pos"]))
  }

})

test_that("SNPinfo stringAsFactor",{
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       verbose = F)

  # simulate SNP
  SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
                                 rep("Chr2", 4),
                                 rep("Chr3", 5)),
                         pos = c(c(1,sample(100, 2)),
                                 c(3,sample(150, 3)),
                                 c(2,sample(200, 4))),
                         SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
                                         1:(3 + 4 + 5)),
                         stringsAsFactors = TRUE)

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  expect_error(expect_identical(SNPs$SNPcoord, SNPcoord))
  expect_is(SNPs$SNPcoord$SNPid, "character")
  expect_is(SNPs$SNPcoord$chr, "character")
  expect_is(SNPs$SNPcoord$pos, "integer")
})






test_that("SNPinfo nSNP method", {
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
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

  expect_equal(SNPs$nSNP(), 3 + 4 + 5)
  expect_named(SNPs$nSNP(1))
  expect_equal(as.numeric(SNPs$nSNP(1)), 3)
  expect_equal(as.numeric(SNPs$nSNP("Chr2")), 4)
  expect_equal(as.numeric(SNPs$nSNP(c("Chr2", "Chr3"))), c(4, 5))


})


test_that("SNPinfo getInfo method", {
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       verbose = F)

  # simulate SNP
  SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
                                 rep("Chr2", 4),
                                 rep("Chr3", 5)),
                         pos = c(c(1,sample(100, 2)),
                                 c(3,sample(150, 3)),
                                 c(2,sample(200, 4))),
                         SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
                                         1:(3 + 4 + 5)),
                         stringsAsFactors = F)

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)


  expect_identical(SNPs$getInfo("SNP01"), SNPs$SNPcoord[1, ])
  expect_equal(SNPs$getInfo(sprintf(fmt = paste0("SNP%0", 2, "i"),
                                    1:(3 + 4 + 5))),
               SNPcoord)
})


test_that("SNPinfo plot method", {
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
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
  p <- SNPs$plot(alpha = 1)

  expect_is(p, "plotly")

})



test_that("specie's \"print\" methods", {
  # create specie
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
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

  expect_output(print(SNPs), "specie: Undefinded")
  expect_output(print(SNPs), "12 Markers on 3 chromosomes :")
  expect_output(print(SNPs), "Chr1 Chr2 Chr3")
  expect_output(print(SNPs), "   3    4    5")
  expect_output(print(SNPs), "SNPcoord:")
  expect_output(print(SNPs), "    chr pos SNPid")
})
