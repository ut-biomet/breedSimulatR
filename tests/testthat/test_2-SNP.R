# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test file for the SNPinfo class


# set.seed(6705) # for reproductible RNG


##### Initialisation functions ####
if (basename(getwd()) == "breedSimulatR") {
  devtools::load_all()
  source("tests/testthat/src/functionsForTests.R")
} else source("src/functionsForTests.R")



#### TESTS ####
test_that("SNPinfo initialization without link map pos", {
  #### Initialisation:
  mySpec <- create_spec()
  ( mySpec <- create_spec() )
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr],
                                replace = FALSE)
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         physPos = physPos,
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
    expect_true(!is.unsorted(SNPs$SNPcoord[SNPs$SNPcoord$chr == chr, "physPos"]))
  }

})

test_that("SNPinfo initialization with link map pos", {
  #### Initialisation:
  mySpec <- create_spec()
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  SNPcoord <- do.call(rbind,lapply(seq(mySpec$nChr),
                           function(chr){
                             physPos <- sample(mySpec$lchr[chr], nMarkers[chr],
                                    replace = FALSE)
                             linkMapPos <- .simulLinkMapPos(physPos,
                                                    mySpec$lchr[chr],
                                                    mySpec$lchrCm[chr],
                                                    b1 = 10, b2 = 10)
                             data.frame(physPos,linkMapPos)
                           }))

  # generate arbitrary SNPid
  SNPcoord$SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord$chr <- rep(mySpec$chrNames, times = nMarkers)

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
    expect_true(!is.unsorted(SNPs$SNPcoord[SNPs$SNPcoord$chr == chr, "physPos"]))
  }
  for (chr in mySpec$chrNames) {
    expect_true(!is.unsorted(SNPs$SNPcoord[SNPs$SNPcoord$chr == chr, "linkMapPos"]))
  }


  # two snpcoord data.frame with same information but in different order should
  # lead to the same SNPinfo object
  isDiff <- FALSE
  while (!isDiff) {
    cols <- sample(colnames(SNPcoord))
    rows <- sample(nrow(SNPcoord))
    if (!all(rows == seq_len(nrow(SNPcoord))) && !all(cols == colnames(SNPcoord))) {
      SNPcoord_2 <- SNPcoord[sample(nrow(SNPcoord)), sample(colnames(SNPcoord))]
      isDiff <- TRUE
    }
  }
  expect_error({SNPs2 <- SNPinfo$new(SNPcoord_2, mySpec)}, NA)
  # browser()
  expect_equal(SNPs2,SNPs)
  expect_identical(SNPs2$SNPcoord, SNPs$SNPcoord)
  expect_identical(SNPs2$ids, SNPs$ids)
  expect_identical(SNPs2$SNPcoordList, SNPs$SNPcoordList)
  expect_identical(capture.output(print(SNPs2)),capture.output(print(SNPs)))
})

test_that("SNPinfo errors", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 5)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))
  # Test differents chromosome names than those in specie
  chrNames <- c("toto", mySpec$chrNames[-1])
  SNPcoord <- data.frame(chr = rep(chrNames, times = nMarkers),
                         physPos = physPos,
                         SNPid = SNPid)
  expect_error({SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)},
               paste('"Chromosomes\'names specified in "SNPcoord"',
                     'do not match those specified in "specie"'))

  SNPcoord <- data.frame(chr = rep(chrNames, times = nMarkers),
                         physPos = physPos,
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
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord_SAF_T <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         physPos = physPos,
                         linkMapPos = NA,
                         SNPid = SNPid,
                         stringsAsFactors = TRUE)
  SNPcoord_SAF_F <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                               physPos = physPos,
                               linkMapPos = NA,
                               SNPid = SNPid,
                               stringsAsFactors = FALSE)

  # create SNPinfo object
  SNPs_SAF_T <- SNPinfo$new(SNPcoord = SNPcoord_SAF_T, specie = mySpec)
  SNPs_SAF_F <- SNPinfo$new(SNPcoord = SNPcoord_SAF_F, specie = mySpec)

  #### Tests:
  expect_error(expect_identical(SNPs_SAF_T$SNPcoord, SNPcoord_SAF_F))
  expect_error(expect_identical(SNPs_SAF_F$SNPcoord, SNPcoord_SAF_F))
  for (col in colnames(SNPcoord_SAF_F)) {
    expect_equal(summary(SNPs_SAF_T$SNPcoord[,col]), summary(SNPcoord_SAF_F[,col]))
  }
  expect_equal(SNPs_SAF_T, SNPs_SAF_F)
})

test_that("SNPinfo nSNP method", {
  #### Initialisation:
  mySpec <- create_spec(nChr = 3)
  nMarkers <- round(mySpec$lchr/10)

  # generate positions
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         physPos = physPos,
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
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         physPos = physPos,
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
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         physPos = physPos,
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
  physPos <- unlist(lapply(seq(mySpec$nChr),
                       function(chr){
                         sample(mySpec$lchr[chr], nMarkers[chr])
                       }))

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarkers)*50, sum(nMarkers)))

  # SNP coordinates data.frame
  SNPcoord <- data.frame(chr = rep(mySpec$chrNames, times = nMarkers),
                         physPos = physPos,
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
