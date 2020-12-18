# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test file for the Specie class

if (interactive()) {
  devtools::load_all()
}


test_that("specie initialization", {

  expect_error({mySpec <- specie$new(nChr = 3,
                                     lchr = c(100, 150, 200),
                                     specName = "Geneticae Exempulus",
                                     ploidy = 2,
                                     mutRate = 10^-8,
                                     recombRate = 10^-7,
                                     chrNames = c("C1", "C2", "C3"),
                                     verbose = F)},
               NA)


  expect_identical(mySpec$specName, "Geneticae Exempulus")
  expect_identical(mySpec$nChr, 3)
  expect_identical(mySpec$ploidy, 2)
  expect_identical(as.numeric(mySpec$lchr), c(100, 150, 200))
  expect_identical(names(mySpec$lchr), c("C1", "C2", "C3"))
  expect_identical(mySpec$mutRate, 10^-8)
  expect_identical(mySpec$chrNames, c("C1", "C2", "C3"))
  expect_identical(mySpec$recombRate, 10^-7)
  expect_output(specie$new(1, 10, verbose = T), paste("A new species has",
                                                      "emerged: Undefinded !"))
  expect_output(specie$new(1, 10, verbose = F), NA)

  expect_error(specie$new(nChr = 3.5,
                          lchr = c(100, 150, 200),
                          verbose = F),
               "nChr must be integer.")
  expect_error(specie$new(nChr = 3,
                          lchr = c(100, 150.2, 200),
                          verbose = F),
               "lchr must be integers.")
})


test_that("specie initialization without optional values", {

  expect_error(specie$new(nChr = 3,
                          lchr = c(100, 150, 200),
                          verbose = F),
               NA)

  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       verbose = F)

  expect_identical(mySpec$specName, "Undefinded")
  expect_identical(mySpec$nChr, 3)
  expect_identical(mySpec$ploidy, 2)
  expect_identical(as.numeric(mySpec$lchr), c(100, 150, 200))
  expect_identical(names(mySpec$lchr), c("Chr1", "Chr2", "Chr3"))
  expect_identical(mySpec$mutRate, NA)
  expect_identical(mySpec$chrNames, c("Chr1", "Chr2", "Chr3"))
  expect_identical(mySpec$recombRate, NA)
})


test_that("specie's \"getChrLength\" methods", {
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       verbose = F)

  expect_is(mySpec$getChrLength(), "numeric")
  expect_named(mySpec$getChrLength())

  expect_identical(as.numeric(mySpec$getChrLength()), c(100, 150, 200))
  expect_identical(as.numeric(mySpec$getChrLength(1)), 100)
  expect_identical(as.numeric(mySpec$getChrLength("Chr2")), 150)

  expect_identical(names(mySpec$getChrLength()), c("Chr1", "Chr2", "Chr3"))
  expect_identical(names(mySpec$getChrLength(3)), "Chr3")
  expect_identical(names(mySpec$getChrLength("Chr1")), "Chr1")
  })


test_that("specie's \"print\" methods", {
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       verbose = F)

  expect_output(print(mySpec), paste("Name: Undefinded\\nNumber of",
                                     "Chromosomes: 3\\nPloidy: 2\\nMutation",
                                     "rate : NA\\nRecombination Rate:",
                                     "NA\\nChromosome length:\\n",
                                     "    chrNames chrLength\\nChr1",
                                     "    Chr1       100\\nChr2",
                                     "    Chr2       150\\nChr3",
                                     "    Chr3       200"))
})
