# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Definition of haplotype class




#### CODE ####

#' R6 Class representing an haplotype
#'
#' @description
#' haplotype object store specific information of individuals haplotype
#'
#'
#' @export
#' @import R6
haplotype <- R6::R6Class(
  "Haplotype",
  public = list(
    #' @field SNPinfo [SNPinfo class] information about the haplotype's SNPs
    #'   (see:\link[breedSimulatR]{SNPinfo})
    SNPinfo = NULL,
    #' @field values [list] list haplotype values
    values = NULL,

    #' @description Create a new Haplotype object.
    #' @param SNPinfo [SNPinfo class] information about the haplotype's SNPs
    #'   (see:\link[breedSimulatR]{SNPinfo})
    #' @param haplo [matrix] named matrix of the genotype for all markers
    #' @return A new `SNPinfo` object.
    #' @examples
    #' # create specie
    #' mySpec <- specie$new(nChr = 3,
    #'                      lchr = c(100, 150, 200),
    #'                      verbose = FALSE)
    #'
    #' # simulate SNP
    #' SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
    #'                                rep("Chr2", 4),
    #'                                rep("Chr3", 5)),
    #'                        pos = c(sample(100, 3),
    #'                                sample(150, 4),
    #'                                sample(200, 5)),
    #'                        SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                                        1:(3 + 4 + 5)))
    #'
    #' # create SNPinfo object
    #' SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)
    #' # simulate haplotype
    #' rawHaplo <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
    #'                    nrow = 2)
    #' colnames(rawHaplo) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                               1:(3 + 4 + 5))
    #' haplo <- haplotype$new(SNPinfo = SNPs,
    #'                        haplo = rawHaplo)
    initialize = function(SNPinfo, haplo) {

      # checks

      # SNPinfo class
      if (class(SNPinfo)[1] != "SNPinfo"){
        stop('class(SNPinfo)[1] != "SNPinfo"\n"SNPinfo" must be a SNPinfo object see: ?SNPinfo')
      }
      # haplo class
      if (class(haplo) != "matrix"){
        stop('class(haplo) != "matrix"\n"haplo" must be a matrix')
      }
      # haplo ploidy
      if (nrow(haplo) != SNPinfo$specie$ploidy) {
        stop('nrow(haplo) != SNPinfo$specie$ploidy\nnrow(haplo) must be equal to the specie ploidy')
      }
      # number of markers
      if (ncol(haplo) != nrow(SNPinfo$SNPcoord)) {
        stop('ncol(haplo) != nrow(SNPinfo$SNPcoord)\nncol(haplo) must be equal to the number of markers in SNPinfo')
      }
      # markers names
      if (is.null(colnames(haplo))) {
        stop('colnames(haplo) is NULL\nhaplo must be a named matrix')
      }
      if (!all(colnames(haplo) %in% SNPinfo$SNPcoord$SNPid)) {
        stop('all(colnames(haplo) %in% SNPinfo$specie$chrNames) is false\ncolnames(haplo) must be the names of the markers')
      }

      self$SNPinfo <- SNPinfo

      self$values <- lapply(SNPinfo$specie$chrNames, function(chr){

        h <- haplo[, SNPinfo$SNPcoord[SNPinfo$SNPcoord$chr == chr, "SNPid"]]
        mode(h) <- "integer"

        rownames(h) <- sprintf(
          fmt = paste0("h%0",
                       floor(log10(SNPinfo$specie$ploidy)) + 1,
                       "i"),
          c(1:SNPinfo$specie$ploidy))
        h
      })
      names(self$values) <- SNPinfo$specie$chrNames

    }),

  active = list(
    #' @field allelDose [numeric] vector of haplotypes encoded in allele dose
    allelDose = function(){
      colSums(do.call(cbind, self$values))
    }

  )
)
