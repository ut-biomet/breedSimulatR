# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of SNPinfo class




#' R6 Class representing a set of SNP markers
#'
#' @description SNPinfo object store specific information of a set of SNP
#' markers
#'
#' @details This class is useful for setting the Haplotype class. (see:
#' \code{\link[breedSimulatR]{haplotype}})
#'
#' @export
#' @importFrom R6 R6Class
SNPinfo <- R6::R6Class(
  "SNPinfo",
  public = list(
    #' @field SNPcoord [data.frame] Coordinate of all SNPs.
    #'
    #' 4 columns:
    #' \itemize{
    #'  \item{\code{chr}:} {Chromosome holding the SNP}
    #'  \item{\code{physPos}:} {SNP physical position on the chromosome}
    #'  \item{\code{linkMapPos}:} {SNP linkage map position on the chromosome}
    #'  \item{\code{SNPid}:} {SNP's IDs}
    #'  }
    SNPcoord = data.frame(
      chr = c(),
      physPos = c(),
      linkMapPos = c(),
      SNPid = c(),
      stringsAsFactors = FALSE
    ),

    #' @field specie [specie class] Specie of the SNPs
    #'   (see:\link[breedSimulatR]{specie})
    specie = NULL,

    #' @field SNPcoordList [list] Named list of dataframes with the coordinate
    #'   of all SNPs. for each chromosomes
    SNPcoordList = list(),

    #' @field ids [list] Named list of the SNP ids for all chromosomes
    ids = list(),

    #' @description Create a new SNPinfo object.
    #' @param SNPcoord [data.frame] Coordinate of all SNPs.
    #'
    #'   3 columns: \itemize{ \item{\code{chr}:} {Chromosome holding the SNP}
    #'   \item{\code{physPos}:} {SNP physical position on the chromosome}
    #'   \item{\code{linkMapPos}:} {SNP linkage map position on the chromosome}
    #'   \item{\code{SNPid}:} {SNP's IDs} }
    #' @param specie [specie class] Specie of the SNPs
    #'   (see:\link[breedSimulatR]{specie})
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
    initialize = function(SNPcoord, specie) {
      # CHECKS:
      # parameters classes
      if (class(specie)[1] != "Specie") {
        stop('"class(specie)" must be "Specie"')
      }
      if (class(SNPcoord) != "data.frame") {
        stop('"class(SNPcoord)" must be "data.frame"')
      }
      # SNPcoord's colnames
      if (ncol(SNPcoord) == 3) {
        if (!all(c("chr", "physPos", "SNPid") %in% colnames(SNPcoord))) {
          stop('"colnames(SNPcoord)" must include "chr", "physPos" and "SNPid", and optionally "linkMapPos"')
        }
        # add linkMapPos column
        SNPcoord$linkMapPos <- NA
      } else if (ncol(SNPcoord) == 4) {
        if (!all(c("chr", "physPos", "SNPid", "linkMapPos") %in% colnames(SNPcoord))) {
          stop('"colnames(SNPcoord)" must include "chr", "physPos" and "SNPid", and optionally "linkMapPos"')
        }
        if (sum(is.na(SNPcoord$linkMapPos)) != 0 &
          sum(is.na(SNPcoord$linkMapPos)) != nrow(SNPcoord)) {
          stop('"SNPcoord$linkMapPos" should either contains only NAs or none')
        }
      }

      # chr names match those in specie
      if (!all(unique(SNPcoord$chr) %in% specie$chrNames)) {
        stop(paste0(
          '"Chromosomes\'names specified in "SNPcoord" do ',
          'not match those specified in "specie"\n',
          "unique(SNPcoord$chr) = ",
          paste0(unique(SNPcoord$chr), collapse = " "), "\n",
          "specie$chrNames = ",
          paste0(specie$chrNames, collapse = " ")
        ))
      }

      # markers are unique
      if (length(SNPcoord$SNPid) != length(unique(SNPcoord$SNPid))) {
        stop('Some markers appear several times in "SNPcoord"')
      }

      # physical position are integers
      if (!all.equal(SNPcoord$physPos, as.integer(SNPcoord$physPos))) {
        stop("markers's physical position should all be integers.")
      }

      # remove factors
      SNPcoord[, vapply(SNPcoord, is.factor, TRUE)] <-
        lapply(SNPcoord[, vapply(SNPcoord, is.factor, TRUE)], as.character)

      # convert to integer
      SNPcoord$physPos <- as.integer(SNPcoord$physPos)

      # sort position in increasing order (needed for the function
      # "findInterval" in individual's generateGametes method)
      SNPcoord <- SNPcoord[order(SNPcoord$SNPid), ] # first order on unique values to be sure to always have the same order.
      SNPcoord <- SNPcoord[order(SNPcoord$physPos), ]

      # reorder columns
      SNPcoord <- SNPcoord[, c("chr", "SNPid", "physPos", "linkMapPos")]
      rownames(SNPcoord) <- SNPcoord$SNPid

      self$SNPcoord <- SNPcoord
      self$specie <- specie
      self$SNPcoordList <- split(
        SNPcoord,
        SNPcoord$chr
      )
      self$SNPcoordList <- self$SNPcoordList[specie$chrNames]

      # check SNP unicity and linkage map position order
      lapply(self$SNPcoordList, function(subSNPcoord) {
        if (!identical(unique(subSNPcoord$physPos), subSNPcoord$physPos)) {
          stop("Some SNPs have the same physical position.")
        }
        # check the linkMapPos are also sorted like the physical position
        if (is.unsorted(subSNPcoord$linkMapPos, na.rm = TRUE, strictly = TRUE)) {
          warning("SNP's position order is not the same when sorted by physical position and by linkage map position.")
        }
        subSNPcoord
      })


      self$ids <- lapply(specie$chrNames, function(chr) {
        SNPcoord[SNPcoord$chr == chr, "SNPid"]
      })
      names(self$ids) <- specie$chrNames
    },

    #' @description
    #' Get the number of SNPs per chromosomes
    #' @param chr [str or numeric] chromosome id
    #' @examples
    #' SNPs$nSNP()
    #' SNPs$nSNP(c("Chr2","Chr3"))
    nSNP = function(chr = NA) {
      if (length(chr) == 1 && is.na(chr)) {
        return(nrow(self$SNPcoord))
      }

      stopifnot((is.character(chr) || is.numeric(chr)))

      if (is.numeric(chr)) {
        chr <- self$specie$chrNames[chr]
      }
      stopifnot(chr %in% self$specie$chrNames)
      vapply(
        chr, function(chr) nrow(self$SNPcoord[self$SNPcoord$chr == chr, ]),
        1
      )
    },
    #' @description
    #' Get information about specific SNPs
    #' @param SNPid [str] SNP ids
    #' @examples
    #' SNPs$getInfo("SNP01")
    #' SNPs$getInfo(c("SNP01", "SNP03"))
    getInfo = function(SNPid) {
      self$SNPcoord[match(SNPid, self$SNPcoord$SNPid), ]
    },

    #' @description Display summary information about the object: specie, number
    #' of SNP, SNP coordinates.
    print = function() {
      cat(paste0(
        "specie: ", self$specie$specName, "\n",
        self$nSNP(), " Markers on ",
        length(unique(self$SNPcoord$chr)), " chromosomes :\n"
      ))
      print(self$nSNP(self$specie$chrNames))
      cat("SNPcoord:\n")
      df <- self$SNPcoord
      df <- df[order(df$SNPid), ]
      print(df)
    },

    #' @description
    #' plot chromosome map using the \pkg{plotly} package
    #' @param alpha transparency see \link[plotly]{plot_ly}
    #'
    #' @examples
    #' SNPs$plot(alpha = 1)
    plot = function(alpha = 0.01) {
      if (requireNamespace("plotly", quietly = TRUE)) {
        ends <- self$specie$lchr

        p <- plotly::plot_ly(
          data = self$SNPcoord,
          x = ~chr,
          y = ~physPos,
          type = "scatter",
          mode = "markers",
          alpha = alpha,
          name = "SNPs",
          hoverinfo = "text",
          text = apply(self$SNPcoord, 1, function(l) {
            paste(names(l), ":", l, collapse = "\n")
          })
        )
        p <- plotly::add_markers(p,
          x = rep(names(ends), 2),
          y = c(ends, rep(0, length(ends))),
          alpha = 1,
          name = "Chromosome's edges",
          hoverinfo = "text",
          text = paste(
            rep(names(ends), 2),
            ": length =",
            rep(ends, 2)
          )
        )
      } else {
        stop("The package 'plotly' is needed to use the method 'plot'.")
      }
    }
  )
)
