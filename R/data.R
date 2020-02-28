# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Data documentation file


#' Example of data
#'
#' Set of data usefull for trying the package.
#'
#' @details This object is a named list containing 3 elements:
#' \itemize{
#'   \item{\code{genotypes} (\code{data.frame}): }{ Genotypic data encoded in allele
#'   doses of 100 fictitious individuals with 3333 SNP markers.
#'   These individuals have 10 chromosomes of length 10^6 bases pairs.}
#'   \item{\code{snpCoord} (\code{data.frame}): }{ Coordinates of the 3333
#'   individuals' markers. This data.frame contains 3 columns: \code{chr},
#'   \code{pos} and \code{SNPid}.}
#'   \item{\code{snpEffects} (\code{numeric}): }{ "true" effects of the 3333
#'   individuals' markers about a fictitious quantitative trait based on an
#'   additive architecture.}
#' }
#'
#' @docType data
#'
#' @source These data are simulated and was generated using the serious game
#'   "PlantBreedGame" available on GitHub:
#'   \url{https://github.com/timflutre/PlantBreedGame}
#'
#'   Flutre, T., Diot, J., and David, J. (2019). PlantBreedGame:
#'   A Serious Game that Puts Students in the Breeder’s Seat. Crop Science.
#'   \href{https://dl.sciencesocieties.org/publications/cs/abstracts/59/4/1374}{DOI 10.2135/cropsci2019.03.0183le}
#'
#' @aliases genotypes snpCoord snpEffects
#'
#' @examples
#' exampleData$genotypes
#' exampleData$snpCoord
#' exampleData$snpEffects
"exampleData"
