# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of phenotyper class



# Iwata example -----------------------------------------------------------

# Assume genetic effects of qtns. The effects are assumed to follow an exponential distribution here.

# lmbd <- 1.0
# is.qtn <- substr(colnames(initPop$genoMat), 1, 3) == "qtn"
# qtn.genotypes <- initPop$genoMat[, is.qtn]
# qtn.eff <- rexp(n.qtn, lmbd) * sample(c(-1, 1), n.qtn, replace = T)
#
#
# # Simulate genotypic values and phenotypic values.
#
# he <- 0.5  # heritability intended
# g <- qtn.genotypes %*% qtn.eff
# var.e <- var(g) * ((1 - he) / he)
# e <- rnorm(length(g), sd = sqrt(var.e))
# y <- g + e
# (var(g) / var(y))  # heritabilty realized
#
# g <- a + d + i
# y <- g + e


# trait -------------------------------------------------------------------

#' R6 class representing a phenotypic trait
#'
#' @description
#' trait object store information about a phenotypic trait
#'
#'
#' @export
#' @import R6
trait <- R6::R6Class(
  "trait",
  lock_objects = FALSE,
  public = list(
    #' @field name [string] Name of the trait
    name = NULL,
    #' @field class "qualitative" or "quantitative"
    class = "quantitative",
    #' @field arch Genetic architecture
    arch = "additive",
    #' @field qtn [character vector] list of the quantitative trait nucleotides
    qtn = NULL,
    #' @field qtnEff [character vector] quantitative trait nucleotides effects
    qtnEff = NULL,
    #' @description Create a new trait object.
    #' @param name [] name description
    #' @param class "quantitative" or "qualitative"
    #' @param arch "arch Genetic architecture
    #' @param qtn qtn [character vector] list of the quantitative trait
    #'   nucleotides
    #' @param qtnEff qtnEff [numeric vector] quantitative trait nucleotides
    #'   effects
    #' @return A new `trait` object.
    #' @details \code{qtn} and \code{qtnEff} must be in the same order.
    #'   \code{qtnEff[n]} must be the QTN effect of the QTN \code{qtn[n]}.
    #' @examples
    #' mySpec <- specie$new(nChr = 10,
    #'                      lchr = 10^6,
    #'                      specName = "Geneticae Exempli",
    #'                      ploidy = 2,
    #'                      recombRate = 3/10^6)
    #' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
    #'                     specie = mySpec)
    #'
    #' myTrait <- trait$new(name = "myTrait",
    #'                      class = "quantitative",
    #'                      qtn = sample(SNPs$SNPcoord$SNPid, 100),
    #'                      qtnEff = rnorm(100, sd = 0.5))
    initialize = function(name = NULL,
                          class = "quantitative",
                          arch  = "additive" ,
                          qtn = NULL,
                          qtnEff = NULL){
      # checks ----
      if (is.null(name)) {
        name <- "Unspecified"
      } else name <- as.character(name)


      if (identical(class, "quantitative")) {
        # checks for "quantitative" traits

        if (identical(arch, "additive")) {
          # checks for "additive" architecture
          if (!is.character(qtn)) {
            stop('"qtn" must be a character vector')
          }
          if (!is.numeric(qtnEff)) {
            stop('"qtnEff" must be a numeric vector')
          }
          if (length(qtn) == 0) {
            stop('"length(qtn)" must be greater than 0')
          }
          if (length(qtn) != length(qtnEff) ) {
            stop('"length(qtn)" must be equal to "length(qtnEff)"')
          }

        } else {
          # not recognized architecture
          stop('"arch" must be "additive", other options not implemented yet')
        }

        # checks for "quantitative" traits
      } else if (identical(arch, "qualitative")) {
        stop('not implemented yet, please try with quantitative trait')
      } else {
        stop('"class" must be "quantitative" or "qualitative"')
      }

      #  trait initialization ----
      self$name <- name
      self$class <- class
      self$arch <- arch
      self$qtn <- qtn
      self$qtnEff <- qtnEff

    },
    #' @description Calculate the genetic values of a population
    #' @param pop [population class] population
    #'   (see:\link[breedSimulatR]{population})
    #' @examples
    #' # create population
    #' example_pop <- createPop(geno = exampleData$genotypes,
    #'                          SNPinfo = SNPs,
    #'                          popName = "Example population")
    #' myTrait$gv(example_pop)
    gv = function(pop){
      # check ----
      if (!identical(class(pop), c("population", "R6"))) {
        stop('"class(pop)" must be a "population"')
      }
      genoMat <- pop$genoMat
      if (!all(self$qtn %in% colnames(genoMat))) {
        stop(paste('qtn are not recognized.',
                   'Please check the QTN are defined in the SNP informations'))
      }

      # gv calculation ----
      if (identical(self$class, "quantitative")) {
        # for "quantitative" traits

        if (identical(self$arch, "additive")) {
          # for "additive" architecture
          gv <- genoMat[, self$qtn] %*% self$qtnEff

        } else {
          # not recognized architecture
          stop('not implemented yet')
        }

        # checks for "quantitative" traits
      } else if (identical(self$arch, "qualitative")) {
        stop('not implemented yet, please try with quantitative trait')
      }

      # out ----
      gv

    },
    #' @description
    #' Display information about the object
    print = function() {
      cat(paste0("trait: ", self$name, "\n",
                 self$class, " trait\n",
                 self$arch, " architecture\n",
                 "number of QTN effects: ", length(self$qtnEff), "\n"))
    }
  )
)
