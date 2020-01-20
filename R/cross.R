# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of individual crossing functions


#' Cross two individuals together
#'
#' @param ind1 parent 1
#' @param ind2 parent 2
#' @param n number of descendants
#' @param names names of the descendants
#' @param verbose print informations
#'
#' @return list of new individuals
#'
#' @examples
makeSingleCross <- function(ind1, ind2, n = 1, names = NULL, verbose = T){

  #  Names
  if (is.null(names)) {
    names <- paste0(ind1$name,"_x_", ind2$name,
                    "_", round(as.numeric(Sys.time())))
    names <- paste0(names, "_", c(1:n))
  }


  gam1 <- ind1$generateGametes(n)
  gam2 <- ind2$generateGametes(n)

  #  Haplotype
  haplo <- mapply(function(g1, g2){
    rbind(g1, g2)
  }, gam1, gam2,
  SIMPLIFY = F)

  newInds <- lapply(c(1:n), function(i){
    individual$new(name = names[[i]],
                   specie = ind1$specie,
                   parent1 = ind1$name,
                   parent2 = ind2$name,
                   haplo = haplotype$new(ind1$haplo$SNPinfo, haplo[[i]]),
                   verbose = F)
  })

  names(newInds) <- names

  newInds
}




#' Proceed to several crosses
#'
#' @param crosses data.frame with crossing instructions: parents names
#' \code{ind1} \code{ind2}, number of descendents \code{n} and names of
#' descendents \code{names}
#' @param pop list of individuals containing the parents
#'
#' @return list of new individuals
#' @export
#'
#' @examples
makeCrosses <- function(crosses, pop){

  newInds <- mapply(makeSingleCross,
                    pop[crosses$ind1], pop[crosses$ind2], crosses$n, crosses$name,
                    USE.NAMES = FALSE )

  newInds
}
