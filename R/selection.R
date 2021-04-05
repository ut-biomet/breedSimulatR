# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of selection functions



#' Selection according to breeding values
#'
#' @param pop (population object) Population of individuals to select
#'   (see:\link[breedSimulatR]{population})
#' @param n (integer) number of individuals to select
#' @param QTNeffects (named numeric vector) effect of the genetic markers
#'
#' @return character vector of the selected individuals' names
#' @export
#'
#' @examples
#' mySpec <- specie$new(specName = "Statisticae exempli",
#'                      nChr = 10,
#'                      lchr = 1e6,
#'                      lchrCm = 100,
#'                      ploidy = 2)
#' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord, specie = mySpec)
#' initPop <- createPop(geno = exampleData$genotypes,
#'                      SNPinfo = SNPs,
#'                      popName = "Initial population")
#' myTrait <- trait$new(name = "myTrait",
#'                      qtn = names(exampleData$snpEffects),
#'                      qtnEff = exampleData$snpEffects)
#' selectBV(pop = initPop,
#'          QTNeffects = myTrait$qtnEff,
#'          n = 10)
selectBV <- function(pop, n, QTNeffects){

  if (is.null(names(QTNeffects))) {
    stop("QTNeffects should be named according to qtns")
  }

  if (!all(names(QTNeffects) %in% colnames(pop$genoMat))) {
    stop("qtn not found in the population. please check that:\n",
         "`all(names(QTNeffects) %in% colnames(pop$genoMat))` is true")
  }

  BV <- as.numeric(pop$genoMat[,names(QTNeffects)] %*% QTNeffects)
  names(pop$inds)[order(BV, decreasing = TRUE)][1:n]
}


#' Selection according to weighted breeding values
#'
#' @param pop (population object) Population of individuals to select
#'   (see:\link[breedSimulatR]{population})
#' @param n (integer) number of individuals to select
#' @param QTNeffects (numeric vector) effect of the genetic markers
#'
#' @return character vector of the selected individuals' names
#' @references Jannink, Jean-Luc. “Dynamics of Long-Term Genomic Selection.”
#'   Genetics Selection Evolution 42, no. 1 (December 2010).
#'   https://doi.org/10.1186/1297-9686-42-35.
#' @export
#'
#' @examples
#' mySpec <- specie$new(specName = "Statisticae exempli",
#'                      nChr = 10,
#'                      lchr = 1e6,
#'                      lchrCm = 100,
#'                      ploidy = 2)
#' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord, specie = mySpec)
#' initPop <- createPop(geno = exampleData$genotypes,
#'                      SNPinfo = SNPs,
#'                      popName = "Initial population")
#' myTrait <- trait$new(name = "myTrait",
#'                      qtn = names(exampleData$snpEffects),
#'                      qtnEff = exampleData$snpEffects)
#' selectWBV(pop = initPop,
#'           QTNeffects = myTrait$qtnEff,
#'           n = 10)
selectWBV <- function(pop, n, QTNeffects){
  if (is.null(names(QTNeffects))) {
    stop("QTNeffects should be named according to qtns")
  }

  if (!all(names(QTNeffects) %in% colnames(pop$genoMat))) {
    stop("qtn not found in the population. please check that:\n",
         "`all(names(QTNeffects) %in% colnames(pop$genoMat))` is true")
  }

  favAllel <- as.numeric(QTNeffects > 0)

  w <- pop$af[names(QTNeffects)]
  w[favAllel == 0] <- 1 - w[favAllel == 0]
  w[w == 0] <- 1 # give weight 1 for fixed alleles
  w <- w ^ (-0.5)

  W_QTNeffects <- QTNeffects * w
  WBV <- as.numeric(pop$genoMat[, names(QTNeffects)] %*% W_QTNeffects)
  names(pop$inds)[order(WBV, decreasing = TRUE)][1:n]
}
