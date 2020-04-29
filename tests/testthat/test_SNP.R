# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test file for the SNPinfo class


# set.seed(6705) # for reproductible RNG


##### Initialisation functions ####
if (interactive()) {
  devtools::load_all()
  source("tests/testthat/src/functionsForTests.R")
} else source("src/functionsForTests.R")



#### TESTS ####
test_that("SNPinfo initialization", {
  #### Initialisation:
  mySpec <- create_spec()
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid)

  #### Tests:
  # create SNPinfo object
  expect_error({SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)},
               NA)

  expect_identical(SNPs$specie, mySpec)
  expect_is(SNPs$SNPcoord, "data.frame")
  expect_is(SNPs$SNPcoordList, "list")
  expect_equal(length(SNPs$SNPcoordList), mySpec$nChr)
  expect_identical(names(SNPs$SNPcoordList), mySpec$chrNames)
  expect_equal(SNPs$nSNP(), sum(nMarkers))

  # check SNPs$SNPcoord is sorted increasingly for all chromosomes (needed for
  # the function "findInterval" in individual's generateGametes method)
  for (chr in mySpec$chrNames) {
    expect_true(!is.unsorted(SNPs$SNPcoord[SNPs$SNPcoord$chr == chr, "pos"]))
  }

})

test_that("SNPinfo errors", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 5)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))
  # Test differents chromosome names than those in specie
  chrNames <- c("toto", mySpec$chrNames[-1])
  SNPcoord <- data.frame(chr = rep(chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid)
  expect_error({SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)},
               paste('"Chromosomes\'names specified in "SNPcoord"',
                     'do not match those specified in "specie"'))

  SNPcoord <- data.frame(chr = rep(chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid,
                         stringsAsFactors = FALSE)
  expect_error({SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)},
               paste('"Chromosomes\'names specified in "SNPcoord"',
                     'do not match those specified in "specie"'))

})




test_that("SNPinfo stringAsFactor",{
  #### Initialisation:
  mySpec <- create_spec()
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord_SAF_T <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid,
                         stringsAsFactors = TRUE)
  SNPcoord_SAF_F <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                               pos = pos,
                               SNPid = SNPid,
                               stringsAsFactors = FALSE)

  # create SNPinfo object
  SNPs_SAF_T <- SNPinfo$new(SNPcoord = SNPcoord_SAF_T, specie = mySpec)
  SNPs_SAF_F <- SNPinfo$new(SNPcoord = SNPcoord_SAF_F, specie = mySpec)

  #### Tests:
  expect_error(expect_identical(SNPs_SAF_T$SNPcoord, SNPcoord_SAF_T))
  expect_equal(summary(SNPs_SAF_T$SNPcoord), summary(SNPcoord_SAF_F))
  expect_equal(SNPs_SAF_T, SNPs_SAF_F)
})

test_that("SNPinfo nSNP method", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 3)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid)

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  #### Tests:
  expect_equal(SNPs$nSNP(), sum(nMarkers))
  expect_named(SNPs$nSNP(1))
  expect_equal(SNPs$nSNP(1), nMarkers[1])
  expect_equal(SNPs$nSNP("Chr2"), nMarkers[2])
  expect_equal(SNPs$nSNP(c("Chr2", "Chr3")), nMarkers[c(2, 3)])
})


test_that("SNPinfo getInfo method", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 3)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid)

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  snp <- sample(SNPid, 1)
  expect_identical(SNPs$getInfo(snp),
                   SNPs$SNPcoord[SNPs$SNPcoord$SNPid == snp, ])
  snps <- sort(sample(SNPid, 10))
  ref <- SNPs$SNPcoord[SNPs$SNPcoord$SNPid %in% snps, ]
  ref <- ref[order(ref$SNPid), ]
  expect_equal(SNPs$getInfo(snps),
               ref)
})


test_that("SNPinfo plot method", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 3)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid)

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)


  expect_error({p <- SNPs$plot(alpha = 1)}, NA)
  expect_is(p, "plotly")
})



test_that("SNPinfo's \"print\" methods", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 20)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  pos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         pos = pos,
                         SNPid = SNPid)

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)

  expect_output(print(SNPs), "specie: Undefinded")
  expect_output(print(SNPs),
                paste(sum(nMarkers),
                      "Markers on",
                      mySpec$nChr, "chromosomes :"))
  expect_output(print(SNPs), "SNPcoord:")
})
