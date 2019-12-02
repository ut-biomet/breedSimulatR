#' R6 Class representing a set of SNP markers
#'
#' @description
#' SNPinfo object store specific information of a set of SNP markers
#'
#' @details
#' This class is usefull for setting the Haplotype class. (see: \code{\link[breedSimulatR]{haplotype}})
#'
#' @export
#' @import R6
SNPinfo <- R6::R6Class(
  "SNPinfo",
  public = list(
    #' @field SNPcoord [data.frame] Coordinate of all SNPs.
    #'
    #' 3 columns:
    #' \itemize{
    #'  \item{\code{chr}:} {Chromosome holding the SNP}
    #'  \item{\code{pos}:} {SNP position on the chromosome}
    #'  \item{\code{SNPid}:} {SNP's IDs}
    #'  }
    SNPcoord = data.frame(chr = c(),
                          pos = c(),
                          SNPid = c(),
                          stringsAsFactors = FALSE),

    #' @field specie [specie class] Specie of the SNPs
    #'   (see:\link[breedSimulatR]{specie})
    specie = NULL,

    #' @field SNPcoordList [list] Named list of dataframes with the coordinate
    #'   of all SNPs. for each chromosomes
    SNPcoordList = list(),

    #' @description
    #' Create a new SNPinfo object.
    #' @param SNPcoord [data.frame] Coordinate of all SNPs.
    #'
    #' 3 columns:
    #' \itemize{
    #'  \item{\code{chr}:} {Chromosome holding the SNP}
    #'  \item{\code{pos}:} {SNP position on the chromosome}
    #'  \item{\code{SNPid}:} {SNP's IDs}
    #'  }
    #' @param specie [specie class] Specie of the SNPs (see:\link[breedSimulatR]{specie})
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
    initialize = function(SNPcoord, specie) {
      stopifnot(all(colnames(SNPcoord) %in% c("chr", "pos", "SNPid")),
                all(levels(SNPcoord$chr) %in% specie$chrNames),
                length(SNPcoord$SNPid) == length(unique(SNPcoord$SNPid)))

      # remove factors
      SNPcoord[, sapply(SNPcoord, is.factor)] <- lapply(SNPcoord[, sapply(SNPcoord, is.factor)], as.character)

      # convert to integer
      SNPcoord$pos <- as.integer(SNPcoord$pos)

      # sort position in increasing order (needed for the function
      # "findInterval" in individual's generateGametes method)
      SNPcoord <- SNPcoord[order(SNPcoord$pos),]


      self$SNPcoord <- SNPcoord
      self$specie <- specie
      self$SNPcoordList <- split(SNPcoord,
                                 SNPcoord$chr)
    },

    #' @description
    #' Get the number of SNPs per chromosomes
    #' @param chr [str or numeric] chromosome id
    #' @examples
    #' SNPs$nSNP()
    #' SNPs$nSNP(c("Chr2","Chr3"))
    nSNP = function(chr=NA) {

      if (length(chr) == 1 && is.na(chr)) return(nrow(self$SNPcoord))

      stopifnot((is.character(chr) || is.numeric(chr)))

      if (is.numeric(chr)) {
        chr <- self$specie$chrNames[chr]
      }
      stopifnot(chr %in% self$specie$chrNames)
      sapply(chr, function(chr) nrow(self$SNPcoord[self$SNPcoord$chr == chr,]))

    },
    #' @description
    #' Get information about specifique SNPs
    #' @param SNPid [str] SNP ids
    #' @examples
    #' SNPs$getInfo("SNP01")
    #' SNPs$getInfo(c("SNP01", "SNP03"))
    getInfo = function(SNPid) {
      self$SNPcoord[match(SNPid, self$SNPcoord$SNPid),]
    },

    #' @description
    #' Display summary information about the object: specie, number of SNP, SNP coordinates.
    print = function() {
      cat(paste0(
        "specie: ", self$specie$specName, "\n",
        self$nSNP(), " Markers on ",
        length(unique(self$SNPcoord$chr)), " chromosomes :\n"
      ))
      print(self$nSNP(self$specie$chrNames))
      cat("SNPcoord:\n")
      print(self$SNPcoord)
    },

    #' @description
    #' plot chromosome map using the \pkg{plotly} package
    #' @param alpha transparency see \link[plotly]{plot_ly}
    #'
    #' @import plotly
    #'
    #' @examples
    #' SNPs$plot(alpha = 1)
    plot = function(alpha = 0.01) {
      ends <- self$specie$lchr

      plotly::plot_ly(data = self$SNPcoord,
                      x = ~chr,
                      y = ~pos,
                      type = "scatter",
                      mode = "markers",
                      alpha = alpha,
                      name = "SNPs",
                      hoverinfo = 'text',
                      text = apply(self$SNPcoord, 1, function(l) {
                        paste(names(l), ":", l, collapse = "\n")
                      })) %>%
        plotly::add_markers(x = rep(names(ends), 2),
                            y = c(ends, rep(0,length(ends))),
                            alpha = 1,
                            name = "Chromosome's edges",
                            hoverinfo = 'text',
                            text = paste(rep(names(ends), 2),
                                         ": length =",
                                         rep(ends, 2)))
    }
  )
)
