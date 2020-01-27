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
  if (is.null(names) || is.na(names)) {
    names <- paste0(ind1$name," x ", ind2$name,
                    "_", as.numeric(Sys.time())*10^5)
    names <- paste0(names, "-", c(1:n))
  }

  if (length(names) != n) {
    if (length(names) == 1) {
      names <- paste0(names, "-", c(1:n))
    } else {
      stop('length(names) should either be equal to one or to n')
    }
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

  # checks
  crosses$ind1 <- as.character(crosses$ind1)
  crosses$ind2 <- as.character(crosses$ind2)
  crosses$names <- as.character(crosses$names)

  # no NA in parnets
  if (any(is.na(crosses$ind1)) || any(is.na(crosses$ind2))) {
    stop('Columns "ind1" and "ind2" should not contain any "NA"')
  }

  # parents in population
  if (any(!crosses$ind1 %in% names(pop$inds))
      || any(!crosses$ind2 %in% names(pop$inds))) {
    notFound <- crosses$ind1[which(! crosses$ind1 %in% names(pop$inds))]
    notFound <- c(notFound,
                  crosses$ind2[which(! crosses$ind2 %in% names(pop$inds))])
    notFound <- paste(notFound, collapse = '" ; "')
    stop(paste0('Parents not found in the population: "', notFound, '"'))
  }


  # number of offspring
  crosses$n[which(is.na(crosses$n))] <- 1
  if (!is.numeric(crosses$n)) {
    stop(paste0('Column n should be numeric values'))
  }
  if (!all(crosses$n == floor(crosses$n))) {
    stop(paste0('Column n should be integer values'))
  }

  newInds <- mapply(makeSingleCross,
                    pop$inds[crosses$ind1],
                    pop$inds[crosses$ind2],
                    crosses$n,
                    crosses$name,
                    USE.NAMES = FALSE )

  newInds <- unlist(newInds)

  # offsprings names do not already exist in population
  if (any(sapply(newInds, function(x){x$name}) %in% names(pop$inds))) {
    nameInPop <- unique(crosses$names[which(crosses$names %in% names(pop$inds))])
    nameInPop <- paste(nameInPop, collapse = '" ; "')
    message(paste0('Offspring names already exist in the population: "', nameInPop, '"'))
  }

  newInds

}
