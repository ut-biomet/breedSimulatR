# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Data documentation file

#' Example of genotypic data.frame
#'
#' A data frame with genotypic data encoded in allele doses: 100 individuals (from specie \link[breedSimulatR]{example_Specie}) and
#' 3333 markers (see: \link[breedSimulatR]{example_SNPs}).
#'
#' @docType data
#' @source These data are simulated and was generated using the serious game
#'   "PlantBreedGame". See: \url{https://github.com/timflutre/PlantBreedGame}
"example_genotypes"

#' Example of Specie object
#'
#' Example of Specie object representing a fake diploid specie "Statisticae
#' exempli" with 10 chromosome of length 10^6.
#'
#' @docType data
#' @source These data are based on simulated data from the serious game
#'   "PlantBreedGame". See: \url{https://github.com/timflutre/PlantBreedGame}
"example_Specie"

#' Example of SNPinfo object
#'
#' Example of SNPinfo object representing SNP information about the
#' \link[breedSimulatR]{example_genotypes} data.
#'
#' @docType data
"example_SNPs"

#' Example of population object
#'
#' Example of population object representing 100 individuals (from specie
#' \link[breedSimulatR]{example_Specie}) the genotypes of these individuals came
#' from \link[breedSimulatR]{example_genotypes} data.
#'
#' @docType data
"example_pop"
