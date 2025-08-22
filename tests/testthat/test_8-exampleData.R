# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Test File -- example data.

##### Initialisation functions ####
if (basename(getwd()) == "breedSimulatR") {
  devtools::load_all()
}



#### TESTS readVCF ####
capture_output({
  test_that("exampleData", {
    expect_equal(names(exampleData), c("snpCoord", "genotypes", "snpEffects"))

    # snp id match
    expect_equal(sort(names(exampleData$snpEffects)), sort(exampleData$snpCoord$SNPid))
    expect_equal(sort(colnames(exampleData$genotypes)), sort(exampleData$snpCoord$SNPid))

    # chr id match new specie's default value
    nChr <- length(unique(exampleData$snpCoord$chr))
    newSpec <- specie$new(lchr = 1e6, lchrCm = 100, nChr = nChr)
    expect_equal(sort(unique(exampleData$snpCoord$chr)), sort(newSpec$chrNames))
  })
})
