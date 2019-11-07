#' Specie class
#'
#' An R6 class defining a specie.
#'
#'@section Class Constructor:
#' \describe{
#'     \item{\code{new(nChr, lchr, specName = "Undefinded", ploidy = NA, mutRate = NA, recombRate = NA,}}{
#'         \itemize{
#'             \item{Generate a new specie}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{nChr}}: [numeric] Number of chromosomes / ploidy
#'                  }
#'                     \item{\bold{\code{lchr}}: [numeric] Length of each chromosomes
#'                  }
#'                     \item{\bold{\code{specName}}: [str] specie's name
#'                  }
#'                     \item{\bold{\code{ploidy}}: [numeric] the number of possible alleles at one locus
#'                  }
#'                     \item{\bold{\code{mutRate}}: [numeric] Mutation rate at each base
#'                  }
#'                     \item{\bold{\code{recombRate}}: [numeric] Mutation rate at each base
#'                  }
#'                     \item{\bold{\code{chrNames}}: [str] Names of the chromosomes
#'                  }
#'                     \item{\bold{\code{verbose}}: [bool]Display info
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{New "specie" object}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#'@section Public Methods:
#' \describe{
#'     \item{\code{getChrLength(chr)}}{
#'         \itemize{
#'             \item{Get the chromosome length}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{chr}}: [str or numeric] chromosome ids
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{Named vector of the chromosomes length}
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#'@section Public Fields:
#' \describe{
#'     \item{\bold{\code{specName}}}{: [str] Specie's name}
#'     \item{\bold{\code{nChr}}}{: [numeric] Number of chromosomes}
#'     \item{\bold{\code{ploidy}}}{: [numeric] Number of possible alleles at one locus}
#'     \item{\bold{\code{lchr}}}{: [numeric] length of all chromosomes}
#'     \item{\bold{\code{mutRate}}}{: [numeric] Mutation rate at each base}
#'     \item{\bold{\code{chrNames}}}{: [str] Names of the chromosomes}
#'     \item{\bold{\code{recombRate}}}{: [numeric] Mutation rate at each base}
#' }
#'
#'
#'@section Special Methods:
#' \describe{
#'     \item{\code{clone(deep = FALSE)}}{
#'         \itemize{
#'             \item{Method for copying an object. See \href{https://adv-r.hadley.nz/r6.html#r6-semantics}{emph{Advanced R}} for the intricacies of R6 reference semantics.}
#'             \item{\bold{Args:}}{
#'                 \itemize{
#'                     \item{\bold{\code{deep}}: logical. Whether to recursively clone nested R6 objects.
#'                  }
#'                 }
#'             }
#'             \item{\bold{Returns:}}{
#'                 \itemize{
#'                     \item{Cloned object of this class.}
#'                 }
#'             }
#'         }
#'     }
#'     \item{\code{print()}}{
#'         \itemize{
#'             \item{Display information about the object}
#'         }
#'     }
#' }
#'
#' @return
#' Specie R6 object
#'
#' @examples
#' mySpec <- specie$new(nChr = 3,
#'                      lchr = c(100, 150, 200),
#'                      specName = "Geneticae Exempulus",
#'                      ploidy = 2,
#'                      mutRate = 10^-8,
#'                      recombRate = 10^-7)
#' print(mySpec)
#'
#' @export
#' @import R6
specie <- R6::R6Class(
  "Specie",
  public = list(
    specName = "Undefinded",
    nChr = NA,
    ploidy = NA,
    lchr = NA,
    mutRate = NA,
    chrNames = NA,
    recombRate = NA,

    initialize = function(nChr,
                          lchr,
                          specName = "Undefinded",
                          ploidy = NA,
                          mutRate = NA,
                          recombRate = NA,
                          chrNames = NA,
                          verbose = T){
      if (!is.numeric(nChr)) stop("nChr must be numeric.")
      if (!is.numeric(lchr)) stop("lchr must be numeric.")
      if (length(lchr) != 1 && length(lchr) != nChr){
        stop(paste("length(lchr) must be equal to 1 (all chr have the same",
                   "size) or equal to nChr."))
      }
      self$specName = specName
      self$nChr = nChr
      self$ploidy = ploidy
      self$mutRate = mutRate
      self$recombRate = recombRate
      if (is.na(chrNames)) {
        self$chrNames <- sprintf(fmt = paste0("Chr%0",
                                              floor(log10(self$nChr)) + 1,
                                              "i"),
                                 1:self$nChr)
      }

      if (length(lchr) == 1) {
        self$lchr = rep(lchr, nChr)
      } else {
        self$lchr = lchr
      }
      names(self$lchr) <- self$chrNames

      if (verbose) cat(paste("A new species has emerged:", self$specName, "!\n\n"))
    },

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