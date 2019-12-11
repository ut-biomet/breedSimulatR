# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Definition of populqtion class




#' R6 Class Representing a Specie
#'
#' @description
#' Specie object store specific information of one specie.
#'
# @details
# Details: Specie object store specific information of one specie.
#'
#' @export
#' @import R6
specie <- R6::R6Class(
  "Specie",
  public = list(
    #' @field specName [str] Specie's name
    specName = "Undefinded",
    #' @field nChr [numeric] Number of chromosomes
    nChr = NA,
    #' @field ploidy [numeric] Number of possible alleles at one locus
    ploidy = NA,
    #' @field lchr [numeric] length of all chromosomes
    lchr = NA,
    #' @field mutRate [numeric] Mutation rate at each base
    mutRate = NA,
    #' @field chrNames [str] Names of the chromosomes
    chrNames = NA,
    #' @field recombRate [numeric] Recombination rate at each base
    recombRate = NA,

    #' @description
    #' Create a new specie object.
    #' @param nChr [numeric] Number of chromosomes
    #' @param lchr [numeric] length of all chromosomes
    #' @param specName [str] Specie's name (optional)
    #' @param ploidy [numeric] Number of possible alleles at one locus (optional)
    #' @param mutRate [numeric] Mutation rate at each base (optional)
    #' @param recombRate [numeric] Recombination rate at each base (optional)
    #' @param chrNames [str] Names of the chromosomes (optional)
    #' @param verbose [bool] Display info (optional)
    #' @return A new `specie` object.
    #' @examples
    #' mySpec <- specie$new(nChr = 3,
    #'                      lchr = c(100, 150, 200),
    #'                      specName = "Geneticae Exempulus",
    #'                      ploidy = 2,
    #'                      mutRate = 10^-8,
    #'                      recombRate = 10^-7)
    #' print(mySpec)
    initialize = function(nChr,
                          lchr,
                          specName = "Undefinded",
                          ploidy = NA,
                          mutRate = NA,
                          recombRate = NA,
                          chrNames = NA,
                          verbose = T){
      if (!is.numeric(nChr)) stop("nChr must be numeric.")
      if (floor(nChr) - nChr != 0) stop("nChr must be integer.")
      if (!is.numeric(lchr)) stop("lchr must be numeric.")
      if (any(floor(lchr) - lchr != 0)) stop("lchr must be integers.")
      if (length(lchr) != 1 && length(lchr) != nChr){
        stop(paste("length(lchr) must be equal to 1 (all chr have the same",
                   "size) or equal to nChr."))
      }
      if (is.na(ploidy)){
        message('"ploidy" was not specify. The ploidy had been set to "2"')
        ploidy <- 2
      }


      self$specName = specName
      self$nChr = nChr
      self$ploidy = ploidy
      self$mutRate = mutRate
      self$recombRate = recombRate
      if (all(is.na(chrNames))) {
        self$chrNames <- sprintf(fmt = paste0("Chr%0",
                                              floor(log10(self$nChr)) + 1,
                                              "i"),
                                 1:self$nChr)
      } else {
        self$chrNames <- chrNames
      }

      if (length(lchr) == 1) {
        self$lchr = rep(lchr, nChr)
      } else {
        self$lchr = lchr
      }
      names(self$lchr) <- self$chrNames

      if (verbose) cat(paste("A new species has emerged:", self$specName, "!\n\n"))
    },

    #' @description
    #' Display information about the object
    print = function() {
      cat(paste0(
        "Name: ", self$specName, "\n",
        "Number of Chromosomes: ", self$nChr, "\n",
        "Ploidy: ", self$ploidy, "\n",
        "Mutation rate : ", self$mutRate, "\n",
        "Recombination Rate: ", self$recombRate, "\n",
        "Chromosome length:\n"
      ))
      print(data.frame(chrNames = self$chrNames,
                       chrLength = self$lchr))
    },

    #' @description
    #' Get the chromosomes length
    #' @param chr [str or numeric] chromosome ids
    #' @examples
    #' mySpec$getChrLength()
    #' mySpec$getChrLength(2)
    #' mySpec$getChrLength("Chr3")
    getChrLength = function(chr = NA) {

      # quick return:
      if (is.na(chr)) {
        return(self$lchr)
      }

      stopifnot((is.character(chr) || is.numeric(chr)))

      if (is.character(chr)) {
        id <- as.numeric(regmatches(chr, gregexpr("[0-9]+", chr)))
      } else id <- chr
      self$lchr[id]

    }
  )
)
