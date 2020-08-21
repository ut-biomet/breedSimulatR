# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of mating functions




#' Random Mating
#'
#' @param inds Vector of the selected individuals' names
#' @param n Total number of cross to perform
#' @param names Names of the new individuals.
#'
#' @return
#' \code{data.frame} of the crossing table
#'
#' @export
#'
#' @examples
#' nCross <- 10
#' inds <- paste0("ind-", seq(1:5))
#' names <- paste0("newInd-", seq(1:nCross))
#' randomMate(inds, nCross, names)
randomMate <- function(inds, n, names) {

  if (length(names) == 1 && n != 1) {
    names <- .charSeq(paste0(names,"-"), seq(n))
  } else if (length(names) != n) {
    stop('"length(names)" must be equal to "1" or to "n"')
  }

  data.frame(
    ind1 = sample(inds, n, replace = TRUE),
    ind2 = sample(inds, n, replace = TRUE),
    n = 1,
    names = names
  )
}
