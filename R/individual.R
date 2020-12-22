# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of individuals class




#' R6 class representing an individual
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
    #' @field haplo [haplotype class] Haplotype of the individual (see:
    #'   \link[breedSimulatR]{haplotype})
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
                          verbose = TRUE){
      self$name <- name
      self$specie <- specie
      self$parent1 <- parent1
      self$parent2 <- parent2
      self$haplo <- haplo
      private$checkHaplo()
      if (verbose) cat(paste("A new ind is borned:", self$name, "!"))
    },

    #' @description Generate Gametes
    #' @param n [float] number of gametes to create (default: 1)
    #' @return list of gametes. A gamete is a named vectors with value 0 or 1.
    #' @examples
    #' myInd$generateGametes()
    #' myInd$generateGametes(2)
    generateGametes = function(n = 1) {
      gametes <- lapply(1:n, function(x) {
        if (any(is.na(self$specie$lchrCm))) {
          stop('specie$lchrCm must be specify in order to generate gemtes')
        }
        gamete <- mapply(function(haplo, SNPcoord) {
          chrName <- SNPcoord$chr[1]
          chrLen <- self$specie$lchr[chrName]
          chrLenCm <- self$specie$lchrCm[chrName]

          # number of recombination locations:
          nRecomb <- rpois(1, chrLenCm/100)

          g <- (runif(1) < 0.5) + 1
          if (nRecomb == 0) {
            gamHaplo <- haplo[g, ]
            return(gamHaplo)
          }



          # draw recombination locations
          if (is.na(SNPcoord$linkMapPos[1])) {
            # if SNP's linkMap positions is unknown
            # draw recombination locations uniformly in physical length range:
            Rpos <- runif(nRecomb, 0, chrLen)
            # get id of the recombination positions
            ids <- (c(0,
                      sort(base::findInterval(Rpos, SNPcoord$physPos)),
                      nrow(SNPcoord))
            )
          } else {
            # draw recombination locations uniformly in chrLenCm range:
            RposCm <- runif(nRecomb, 0, chrLenCm)
            # get id of the recombination positions
            ids <- (c(0,
                      sort(base::findInterval(RposCm, SNPcoord$linkMapPos)),
                      nrow(SNPcoord))
            )
          }

          gamHaplo <- integer(ncol(haplo))
          for (i in seq_len(length(ids) - 1)) {
            if (ids[i] == ids[i + 1]) {
              # several recombinations between SNPs
              g <- -g + 3
              next
            }
            gamHaplo[(ids[i] + 1):ids[i + 1]] <-
              haplo[g, (ids[i] + 1):ids[i + 1]]
            g <- -g + 3
          }

          names(gamHaplo) <- colnames(haplo)
          gamHaplo
        },
        self$haplo$values,
        self$haplo$SNPinfo$SNPcoordList,
        SIMPLIFY = FALSE)

        gamete <- unlist(gamete, use.names = FALSE)
        names(gamete) <- names(self$haplo$allelDose)
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
