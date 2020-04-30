# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of individual crossing functions


#' Cross two individuals together
#'
#' @param ind1 parent 1
#' @param ind2 parent 2
#' @param names names of the descendants
#' @param n number of descendants
#' @param verbose print informations
#'
#' @return list of new individuals
#' @export
#'
#' @examples
#' # create specie:
#' mySpec <- specie$new(nChr = 3,
#'                      lchr = c(100, 150, 200),
#'                      recombRate = 0.006,
#'                      verbose = FALSE)
#'
#' # simulate SNP:
#' SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
#'                                rep("Chr2", 4),
#'                                rep("Chr3", 5)),
#'                        pos = c(sample(100, 3),
#'                                sample(150, 4),
#'                                sample(200, 5)),
#'                        SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                                        1:(3 + 4 + 5)))
#'
#' # create SNPinfo object:
#' SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)
#'
#'
#' # simulate haplotype:
#' rawHaplo1 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
#'                     nrow = 2)
#' colnames(rawHaplo1) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                               1:(3 + 4 + 5))
#' haplo1 <- haplotype$new(SNPinfo = SNPs,
#'                        haplo = rawHaplo1)
#' rawHaplo2 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
#'                     nrow = 2)
#' colnames(rawHaplo2) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                               1:(3 + 4 + 5))
#' haplo2 <- haplotype$new(SNPinfo = SNPs,
#'                        haplo = rawHaplo2)
#'
#'
#' # create individuals:
#' myInd1 <-  individual$new(name = "Ind 1",
#'                          specie = mySpec,
#'                          parent1 = "OkaaSan",
#'                          parent2 = "OtouSan",
#'                          haplo = haplo1,
#'                          verbose = FALSE)
#' myInd2 <-  individual$new(name = "Ind 2",
#'                          specie = mySpec,
#'                          parent1 = "OkaaSan",
#'                          parent2 = "OtouSan",
#'                          haplo = haplo2,
#'                          verbose = FALSE)
#' offspring <- makeSingleCross(myInd1, myInd2, names = "off 1")
#' offspring
makeSingleCross <- function(ind1, ind2, names, n = 1, verbose = TRUE){

  #  Names
  if (is.null(names) || is.na(names)) {
    stop('"names" should be provided')
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
  SIMPLIFY = FALSE)

  newInds <- lapply(c(1:n), function(i){
    individual$new(name = names[[i]],
                   specie = ind1$specie,
                   parent1 = ind1$name,
                   parent2 = ind2$name,
                   haplo = haplotype$new(ind1$haplo$SNPinfo, haplo[[i]]),
                   verbose = FALSE)
  })

  names(newInds) <- names

  newInds
}




#' Proceed several crosses
#'
#' @param crosses data.frame with crossing instructions: parents names
#' \code{ind1} \code{ind2}, number of descendant \code{n} and names of
#' descendant \code{names}
#' @param pop list of individuals containing the parents
#'
#' @return list of new individuals
#' @export
#'
#' @examples
#' # create specie:
#' mySpec <- specie$new(nChr = 3,
#'                      lchr = c(100, 150, 200),
#'                      recombRate = 0.006,
#'                      verbose = FALSE)
#'
#' # simulate SNP:
#' SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
#'                                rep("Chr2", 4),
#'                                rep("Chr3", 5)),
#'                        pos = c(sample(100, 3),
#'                                sample(150, 4),
#'                                sample(200, 5)),
#'                        SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                                        1:(3 + 4 + 5)))
#'
#' # create SNPinfo object:
#' SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)
#'
#'
#' # simulate haplotype:
#' rawHaplo1 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
#'                     nrow = 2)
#' colnames(rawHaplo1) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                               1:(3 + 4 + 5))
#' haplo1 <- haplotype$new(SNPinfo = SNPs,
#'                        haplo = rawHaplo1)
#' rawHaplo2 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
#'                     nrow = 2)
#' colnames(rawHaplo2) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                               1:(3 + 4 + 5))
#' haplo2 <- haplotype$new(SNPinfo = SNPs,
#'                        haplo = rawHaplo2)
#' rawHaplo3 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
#'                     nrow = 2)
#' colnames(rawHaplo3) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
#'                               1:(3 + 4 + 5))
#' haplo3 <- haplotype$new(SNPinfo = SNPs,
#'                        haplo = rawHaplo3)
#'
#'
#' # create individuals:
#' myInd1 <-  individual$new(name = "Ind 1",
#'                          specie = mySpec,
#'                          parent1 = "OkaaSan",
#'                          parent2 = "OtouSan",
#'                          haplo = haplo1,
#'                          verbose = FALSE)
#' myInd2 <-  individual$new(name = "Ind 2",
#'                          specie = mySpec,
#'                          parent1 = "OkaaSan",
#'                          parent2 = "OtouSan",
#'                          haplo = haplo2,
#'                          verbose = FALSE)
#' myInd3 <-  individual$new(name = "Ind 3",
#'                          specie = mySpec,
#'                          parent1 = "OkaaSan",
#'                          parent2 = "OtouSan",
#'                          haplo = haplo3,
#'                          verbose = FALSE)
#' myPop <- population$new(name = "My Population 1",
#'                         inds = list(myInd1, myInd2, myInd3),
#'                         verbose = FALSE)
#'
#' crossToDo <- data.frame(ind1 = c("Ind 1", "Ind 1", "Ind 2"),
#'                          ind2 = c("Ind 2", "Ind 3", "Ind 3"),
#'                          n = 1,
#'                          names = c("Off 1-2", "Off 1-3", "Off 2-3"))
#'
#' makeCrosses(crossToDo, myPop)
makeCrosses <- function(crosses, pop){

  # checks
  if (!all(colnames(crosses)%in%c("ind1", "ind2", "n", "names"))) {
    stop('colnames(crosses) must be "ind1", "ind2", "n", "names".')
  }
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
  # names
  if (any(is.na(crosses$name))) {
    noNames <- crosses[is.na(crosses$name),]
    noNames$name <- paste(noNames$ind1, "x", noNames$ind2, "-",
                          seq_len(nrow(noNames)), "-")
    crosses[is.na(crosses$name),"names"] <- noNames$name
  }

  newInds <- mapply(makeSingleCross,
                    pop$inds[crosses$ind1],
                    pop$inds[crosses$ind2],
                    crosses$name,
                    crosses$n,
                    USE.NAMES = FALSE )

  newInds <- unlist(newInds)

  # offsprings names do not already exist in population
  if (any(vapply(newInds, function(x){x$name}, "character")%in%names(pop$inds))){
    nameInPop <-
      unique(crosses$names[which(crosses$names %in% names(pop$inds))])
    nameInPop <- paste(nameInPop, collapse = '" ; "')
    message(paste0(
      'Offspring names already exist in the population: "',
      nameInPop,
      '"'
    ))
  }

  newInds

}
