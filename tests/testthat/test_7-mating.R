# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Test File -- mating functions.


# set.seed(9552) # for reproductible RNG

if (interactive()) {
  devtools::load_all()
}

validCrossTable <- function(crossTable){
  expect_is(crossTable, "data.frame")
  expect_equal(colnames(crossTable), c("ind1", "ind2", "n", "names"))
  expect_true(!any(is.na(crossTable$ind1)))
  expect_true(!any(is.na(crossTable$ind2)))
  expect_true(all(crossTable$n %% 1 == 0))
  expect_equal(unique(crossTable$names), crossTable$names)
}

#### TESTS randomMate ####
test_that("randomMate", {

  nCross <- sample(1:100, 1)
  inds <- paste0("ind-", seq(1:5))
  names <- paste0("newInd-", seq(1:nCross))

  expect_error({crossTable <- randomMate(inds, nCross, names)}, NA)
  validCrossTable(crossTable)
  expect_true(all(crossTable$ind1 %in% inds))
  expect_true(all(crossTable$ind2 %in% inds))
  expect_equal(as.character(crossTable$names), names)

  names <- "newInd"
  expect_error({crossTable <- randomMate(inds, nCross, names)}, NA)
  validCrossTable(crossTable)

})
