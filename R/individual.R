# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Definition of individuals class




#' An R6 class representing an individual
#'
#' @description
#' individual object store specific information about one individual
#'
#'
#' @export
#' @import R6
individual <- R6::R6Class(
  "individual",
  public = list(
    #' @field name [string] Name of the individual
    name = NULL,
    #' @field specie [specie class] Specie of the SNPs
    #'   (see:\link[breedSimulatR]{specie})
    specie = NULL,
    #' @field parent1 [string] Name of the individual's parent
    parent1 = NULL,
    #' @field parent2 [string] Name of the individual's parent
    parent2 = NULL,
    #' @field haplo [haplotype class] Haplotype of the individual (see: \link[breedSimulatR]{haplotype})
    haplo = NULL,

    #' @description Create a new individual object.
    #' @param name [string] name of the individual
    #' @param specie [specie class] Specie of the SNPs
    #'   (see:\link[breedSimulatR]{specie})
    #' @param parent1 [string] Name of the individual's parent
    #' @param parent2 [string] Name of the individual's parent
    #' @param haplo [haplotype class] Haplotype of the individual (see:
    #'   \link[breedSimulatR]{haplotype})
    #' @param verbose [boolean] display information
    #' @return A new `individual` object.
    #' @examples
    #' # create specie
    #' mySpec <- specie$new(nChr = 3,
    #'                      lchr = c(100, 150, 200),
    #'                      recombRate = 0.006,
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
    #' rawHaplo <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE), nrow = 2)
    #' colnames(rawHaplo) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                               1:(3 + 4 + 5))
    #' haplo <- haplotype$new(SNPinfo = SNPs,
    #'                        haplo = rawHaplo)
    #' myInd <-  individual$new(name = "Ind 1",
    #'                          specie = mySpec,
    #'                          parent1 = "OkaaSan",
    #'                          parent2 = "OtouSan",
    #'                          haplo = haplo,
    #'                          verbose = FALSE)
    initialize = function(name = "Unnamed",
                          specie = specie$new(),
                          parent1 = NA,
                          parent2 = NA,
                          haplo = NA,
                          verbose = T){
      self$name = name
      self$specie = specie
      self$parent1 = parent1
      self$parent2 = parent2
      self$haplo = haplo
      private$checkHaplo()
      if (verbose) cat(paste("A new ind is borned:", self$name, "!"))
    },

    #' @description Get the number of SNPs per chromosomes
    #' @param n [float] number of gemetes to create (default: 1)
    #' @return list of gametes. A gamete is a named vectors with value 0 or 1.
    #' @examples
    #' myInd$generateGametes()
    #' myInd$generateGametes(2)
    generateGametes = function(n = 1) {
      gametes <- lapply(1:n, function(x) {
        if (is.na(self$specie$recombRate)) {
          stop('specie$recombRate must be specify in order to generate gemtes')
        }
        gamete <- mapply(function(haplo, SNPcoord) {
          chrName <- SNPcoord$chr[1]
          chrLen <- self$specie$lchr[chrName]

          # number of recombination events:
          nRecomb <- rbinom(1, chrLen, self$specie$recombRate)

          # recombination place (optimized)
          if (nRecomb <= chrLen/2 && chrLen > 1e7) {
            Rpos <- .Internal(sample2(chrLen, nRecomb))
          } else {
            Rpos <- .Internal(sample(chrLen, nRecomb, F, NULL))
          }

          gamHaplo <- integer(ncol(haplo))
          # split SNP beetween two chromosome
          if (length(Rpos) == 0) {
            g <- (runif(1) < 0.5) + 1
            gamHaplo <- haplo[g, ]
          } else {
            ids <- unique(c(0,
                            sort(findInterval(Rpos, SNPcoord$pos)),
                            length(SNPcoord$pos))
            )
            g <- (runif(1) < 0.5) + 1
            for (i in seq_len(length(ids) - 1)) {
              gamHaplo[(ids[[i]] + 1):ids[[i + 1]]] <-
                haplo[g, (ids[[i]] + 1):ids[[i + 1]]]
              g <- -g + 3
            }
          }
          names(gamHaplo) <- colnames(haplo)
          gamHaplo
        },
        self$haplo$values,
        self$haplo$SNPinfo$SNPcoordList,
        SIMPLIFY = F)

        gamete <- do.call(c, gamete)
        names(gamete) <- sub(".*(?=\\.).", "", names(gamete), perl = TRUE)
        gamete

      })

      gametes
    }
  ),
  private = list(
    # @description
    # Check the number of chromosomes between `haplo` and `specie` is the
    # same
    # @return boolean
    checkHaplo = function(){
      if (length(self$haplo$values) != self$specie$nChr) {
        stop(paste0("Number of chromosomes in haplo is; ",
                    length(self$haplo$values), " but must be equal to: ",
                    self$specie$nChr, " (number of chr of ",
                    self$specie$name, ")"))
      }
    }
  )
)
