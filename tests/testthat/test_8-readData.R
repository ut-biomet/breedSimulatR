# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Test File -- readData functions.

##### Initialisation functions ####
if (basename(getwd()) == "breedSimulatR") {
  devtools::load_all()
  pathPref <- "tests/testthat/"
} else {
  pathPref <- ""
}

phasedVcfFile <- paste0(pathPref, "data/riceDiv_44kSNP_100colls_beagleOut.vcf.gz")
# phasedVcfFile <- paste0(pathPref, "data/breedSim.vcf.gz")
# unPhasedVcfFile <- paste0(pathPref, "data")
# mixedVcfFile <- paste0(pathPref, "data")

#### TESTS readVCF ####
capture_output({
  test_that("readVCF", {
    expect_error({
      dta <- readVCF(file = phasedVcfFile, verbose = F)
    }, NA)
    expect_is(dta, "list")
    expect_identical(names(dta), c("pop", "snps", "specie"))

    # use gaston as reference
    x <- gaston::read.vcf(phasedVcfFile, verbose = FALSE)
    expect_equal(sort(dta$snps$SNPcoord$SNPid), sort(x@snps$id)) 
    expect_equal(sort(names(dta$pop$inds)), sort(x@ped$id)) 
    expect_equal(dta$pop$genoMat, gaston::as.matrix(x))
  })
})
