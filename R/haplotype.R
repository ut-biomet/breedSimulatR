# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of haplotype class




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
    #' @field SNPinfo [SNPinfo class] information about the haplotype of the
    #'   SNPs (see:\link[breedSimulatR]{SNPinfo})
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
    #'                      lchrCm = 100,
    #'                      verbose = FALSE)
    #'
    #' # simulate SNP
    #' SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
    #'                                rep("Chr2", 4),
    #'                                rep("Chr3", 5)),
    #'                        physPos = c(sample(100, 3),
    #'                                sample(150, 4),
    #'                                sample(200, 5)),
    #'                        linkMapPos = NA,
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
      if (class(SNPinfo)[1] != "SNPinfo") {
        stop(paste('class(SNPinfo)[1] != "SNPinfo"\n"SNPinfo" must be a',
             'SNPinfo object see: ?SNPinfo'))
      }
      # haplo class
      if (!is(haplo,"matrix")) {
        stop('is(haplo,"matrix") is FALSE\n"haplo" must be a matrix')
      }
      # haplo ploidy
      if (nrow(haplo) != SNPinfo$specie$ploidy) {
        stop(
          paste('nrow(haplo) != SNPinfo$specie$ploidy\nnrow(haplo) must be',
             'equal to the specie ploidy'))
      }
      # number of markers
      if (ncol(haplo) != nrow(SNPinfo$SNPcoord)) {
        stop(paste(
          'ncol(haplo) != nrow(SNPinfo$SNPcoord)\nncol(haplo) must be equal to',
          'the number of markers in SNPinfo'))
      }
      # markers names
      if (is.null(colnames(haplo))) {
        stop('colnames(haplo) is NULL\nhaplo must be a named matrix')
      }
      if (!all(colnames(haplo) %in% SNPinfo$SNPcoord$SNPid)) {
        stop(paste(
          'all(colnames(haplo) %in% SNPinfo$specie$chrNames)',
          'is false\ncolnames(haplo) must be the names of the markers'
        ))
      }

      self$SNPinfo <- SNPinfo

      self$values <- list()
      for (chr in SNPinfo$specie$chrNames) {
        h <- haplo[, SNPinfo$ids[[chr]]]
        mode(h) <- "integer"
        rownames(h) <- .charSeq("h", c(1:SNPinfo$specie$ploidy))
        self$values[[chr]] <- h
      }

    }),

  active = list(
    #' @field allelDose [numeric] vector of haplotypes encoded in allele dose
    allelDose = function(){
      colSums(do.call(cbind, self$values))
    }

  )
)
