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
                                     lchrCm = c(101, 151, 201),
                                     specName = "Geneticae Exempulus",
                                     chrNames = c("C1", "C2", "C3"),
                                     verbose = F)},
               NA)


  expect_identical(mySpec$specName, "Geneticae Exempulus")
  expect_identical(mySpec$nChr, 3)
  expect_identical(mySpec$ploidy, 2)
  expect_identical(as.numeric(mySpec$lchr), c(100, 150, 200))
  expect_identical(names(mySpec$lchr), c("C1", "C2", "C3"))
  expect_identical(as.numeric(mySpec$lchrCm), c(101, 151, 201))
  expect_identical(names(mySpec$lchrCm), c("C1", "C2", "C3"))
  expect_identical(mySpec$chrNames, c("C1", "C2", "C3"))
  expect_output(specie$new(1, 10, 100, verbose = T), paste("A new species has",
                                                      "emerged: Undefinded !"))
  expect_output(specie$new(1, 10, 100, verbose = F), NA)

  expect_error(specie$new(nChr = 3.5,
                          lchr = c(100, 150, 200),
                          lchrCm = c(101, 151, 201),
                          verbose = F),
               "nChr must be integer.")
  expect_error(specie$new(nChr = 3,
                          lchr = c(100, 150.2, 200),
                          lchrCm = c(101, 151, 201),
                          verbose = F),
               "lchr must be integers.")
})


test_that("specie initialization without optional values", {

  expect_error(specie$new(nChr = 3,
                          lchr = c(100, 150, 200),
                          lchrCm = c(101, 151, 201),
                          verbose = F),
               NA)

  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       lchrCm = c(101, 151, 201),
                       verbose = F)

  expect_identical(mySpec$specName, "Undefinded")
  expect_identical(mySpec$nChr, 3)
  expect_identical(mySpec$ploidy, 2)
  expect_identical(as.numeric(mySpec$lchr), c(100, 150, 200))
  expect_identical(names(mySpec$lchr), c("Chr1", "Chr2", "Chr3"))
  expect_identical(as.numeric(mySpec$lchrCm), c(101, 151, 201))
  expect_identical(names(mySpec$lchrCm), c("Chr1", "Chr2", "Chr3"))
  expect_identical(mySpec$chrNames, c("Chr1", "Chr2", "Chr3"))
})


test_that("specie's \"getChrLength\" methods", {
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       lchrCm = c(101, 151, 201),
                       verbose = F)

  expect_is(mySpec$getChrLength(), "numeric")
  expect_named(mySpec$getChrLength())

  expect_identical(as.numeric(mySpec$getChrLength()), c(100, 150, 200))
  expect_identical(as.numeric(mySpec$getChrLength(1)), 100)
  expect_identical(as.numeric(mySpec$getChrLength("Chr2")), 150)

  expect_identical(names(mySpec$getChrLength()), c("Chr1", "Chr2", "Chr3"))
  expect_identical(names(mySpec$getChrLength(3)), "Chr3")
  expect_identical(names(mySpec$getChrLength("Chr1")), "Chr1")

  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       lchrCm = c(101, 151, 201),
                       verbose = F,
                       chrNames = c("chr1","chr2","X"))
  expect_identical(as.numeric(mySpec$getChrLength("X")), 200)



  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       lchrCm = c(101, 151, 201),
                       verbose = F)

  expect_is(mySpec$getChrLengthCm(), "numeric")
  expect_named(mySpec$getChrLengthCm())

  expect_identical(as.numeric(mySpec$getChrLengthCm()), c(101, 151, 201))
  expect_identical(as.numeric(mySpec$getChrLengthCm(1)), 101)
  expect_identical(as.numeric(mySpec$getChrLengthCm("Chr2")), 151)

  expect_identical(names(mySpec$getChrLengthCm()), c("Chr1", "Chr2", "Chr3"))
  expect_identical(names(mySpec$getChrLengthCm(3)), "Chr3")
  expect_identical(names(mySpec$getChrLengthCm("Chr1")), "Chr1")

  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       lchrCm = c(101, 151, 201),
                       verbose = F,
                       chrNames = c("chr1","chr2","X"))
  expect_identical(as.numeric(mySpec$getChrLengthCm("X")), 201)


  })




test_that("specie's \"print\" methods", {
  mySpec <- specie$new(nChr = 3,
                       lchr = c(100, 150, 200),
                       lchrCm = c(101, 151, 201),
                       verbose = F)

  expect_output(print(mySpec), paste0("Name: Undefinded\\n",
                                      "Number of Chromosomes: 3\\n",
                                      "Ploidy: 2\\n",
                                      "Chromosome length:\\n",
                                      "     chrNames chrLength chrLengthCm\\n",
                                      "Chr1     Chr1       100         101\\n",
                                      "Chr2     Chr2       150         151\\n",
                                      "Chr3     Chr3       200         201"))
})
