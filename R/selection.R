# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of selection functions



#' Selection according to breeding values
#'
#' @param pop (population object) Population of individuals to select
#'   (see:\link[breedSimulatR]{population})
#' @param SNPeffects (numeric vector) effect of the genetic markers
#' @param n (integer) number of individuals to select
#'
#' @return character vector of the selected individulas' names
#' @export
#'
#' @examples
#' mySpec <- specie$new(specName = "Statisticae exempli",
#'                      nChr = 10,
#'                      lchr = 1e6,
#'                      ploidy = 2,
#'                      recombRate = 3/1e6)
#' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord, specie = specie_statEx)
#' initPop <- createPop(geno = exampleData$genotypes,
#'                      SNPinfo = SNPs,
#'                      popName = "Initial population")
#' selectBV(pop = initPop,
#'          SNPeffects = exampleData$snpEffects,
#'          n = 10)
selectBV <- function(pop, SNPeffects, n){
  SNPeffects <- SNPeffects[colnames(pop$genoMat)]
  BV <- as.numeric(pop$genoMat %*% SNPeffects)
  names(pop$inds)[order(BV)][1:n]
}
