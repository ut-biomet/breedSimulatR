# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of phenotyper class


# trait ----

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
    ## Public fields ####
    #' @field name [string] Name of the trait
    name = NULL,
    #' @field class "qualitative" or "quantitative" trait ("qualitative" not
    #' implemented yet)
    class = "quantitative",
    #' @field qtn [character vector] names of the causal quantitative trait
    #' nucleotides (\code{length(qtn) == length(qtnEff)} must be \code{true})
    qtn = NULL,
    #' @field qtnEff [numeric vector] quantitative trait nucleotides effects
    qtnEff = NULL,


    ## Public methods ####
    #' @description Create a new trait object.
    #' @param name [character] name of the trait
    #' @param class "quantitative" or "qualitative"
    #' @param qtn [character vector] list of the quantitative trait
    #'   nucleotides names implied in the trait
    #' @param qtnEff [numeric vector] quantitative trait nucleotides effects
    #' (see details for more information).
    #'
    #' @return A new `trait` object.
    #'
    #' @details \code{qtn} and \code{qtnEff} must be in the same order.
    #'   \code{qtnEff[n]} must be the QTN effect of the QTN \code{qtn[n]}.
    #'
    #'
    #'
    #' @examples
    #' mySpec <- specie$new(nChr = 10,
    #'                      lchr = 10^6,
    #'                      lchrCm = 100,
    #'                      specName = "Geneticae Exempli")
    #' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
    #'                     specie = mySpec)
    #'
    #' myTrait <- trait$new(name = "myTrait",
    #'                      qtn = sample(SNPs$SNPcoord$SNPid, 100),
    #'                      qtnEff = rnorm(100, sd = 0.5))
    initialize = function(name = NULL,
                          class = "quantitative",
                          qtn = NULL,
                          qtnEff = NULL){
      # checks ----
      if (is.null(name)) {
        name <- "Unspecified"
      } else name <- as.character(name)


      if (identical(class, "quantitative")) {
        # checks for "quantitative" traits

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
        if (!is.null(names(qtnEff))) {
          if (any(sort(names(qtnEff)) != sort(qtn))) {
            stop('mismatch between "qtn" and "names(qtnEff)"')
          } else {
            qtnEff <- qtnEff[qtn] # reorder qtnEff
          }
        } else {
          names(qtnEff) <- qtn
        }



        # checks for "quantitative" traits
      } else if (identical(self$class, "qualitative")) {
        stop('not implemented yet, please try with quantitative trait')
      } else {
        stop('"class" must be "quantitative" or "qualitative"')
      }

      #  trait initialization ----
      self$name <- name
      self$class <- class
      self$qtn <- qtn
      self$qtnEff <- qtnEff

    },
    #' @description Calculate the genetic values of a population
    #' @param pop [population class] population
    #'   (see: \link[breedSimulatR]{population})
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
        gv <- genoMat[, self$qtn] %*% self$qtnEff


        # checks for "quantitative" traits
      } else if (identical(self$class, "qualitative")) {
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
                 "number of QTN: ", length(self$qtnEff), "\n"))
    }
  )
)




# phenotyper ----

# TODO find better name than "phenotyper"

#' R6 class representing a phenotyper
#'
#' @description
#' phenotyper object is used to phenotype populations
#'
#' @export
#' @import R6
phenotyper <- R6::R6Class(
  "phenotyper",
  lock_objects = FALSE,
  public = list(
    ## Public Fields ####
    #' @field name [character] name of the phenotyper
    name = "Pheno lab",
    #' @field plotCost [numeric] cost for phenotyping one plot
    plotCost = 1,
    #' @field traits [list] of phenotyped traits
    traits = NULL,
    #' @field mu [numeric] vector of the mean for each trait
    mu = 0,
    #' @field ve [numeric] vector of the environmental variance sigma^2
    #' for each trait
    ve = 1,

    ## Public Methods ####
    #' @description Create a new phenotyper object.
    #' @param name [character] name of the phenotyper
    #' @param traits [trait or list] of phenotyped traits
    #' @param plotCost [numeric] cost for phenotyping one plot
    #' @param mu [numeric] vector of the mean for each trait
    #' @param ve [numeric] (if \code{he} not provided) vector of the
    #' environmental variance sigma^2 for each trait
    #' @param he [numeric] (if \code{ve} not provided) vector of the target
    #' heredity for each trait
    #' @param pop [population] (if \code{ve} not provided) population used for
    #' the target heredity calculation (see:\link[breedSimulatR]{population})
    #' @return A new `phenotyper` object.
    #' @examples
    #' mySpec <- specie$new(nChr = 10,
    #'                     lchr = 10^6,
    #'                     lchrCm = 100,
    #'                     specName = "Geneticae Exempli")
    #' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
    #'                    specie = mySpec)
    #' example_pop <- createPop(geno = exampleData$genotypes,
    #'                         SNPinfo = SNPs,
    #'                         popName = "Example population")
    #' myTrait1 <- trait$new(name = "Trait1",
    #'                      qtn = sample(SNPs$SNPcoord$SNPid, 100),
    #'                      qtnEff = rnorm(100, sd = 0.3))
    #' myTrait2 <- trait$new(name = "Trait2",
    #'                      qtn = sample(SNPs$SNPcoord$SNPid, 100),
    #'                      qtnEff = rnorm(100, sd = 0.4))
    #'
    #' phenoLab1 <- phenotyper$new(name = "My phenoLab",
    #'                             traits = list(myTrait1, myTrait2),
    #'                             plotCost = 1,
    #'                             mu = c(100, 50),
    #'                             ve = NULL,
    #'                             he = c(0.4, 0.55),
    #'                             pop = example_pop)
    #'
    #' phenoLab2 <- phenotyper$new(name = "My phenoLab",
    #'                            traits = list(myTrait1, myTrait2),
    #'                            plotCost = 2,
    #'                            mu = c(100, 50),
    #'                            ve = c(9, 7),
    #'                            he = NULL,
    #'                            pop = NULL)
    #'

    initialize = function(name,
                          traits,
                          plotCost = 1,
                          mu = 0,
                          ve = NULL,
                          he = NULL,
                          pop = NULL){

      if (!is.character(name)) {
        stop("name should be a character string.")
      }

      if (!is.list(traits) && class(traits)[1] != "trait") {
        stop("traits should be list of traits or a single trait object")
      } else if (is.list(traits)) {
        if (length(traits) < 1) {
          stop("The list of traits should not be empty")
        }
        if (!all(vapply(traits, function(x) {class(x)[1]},
                        FUN.VALUE = "trait") == "trait")) {
          stop("All elements of the trait list should be traits")
        }
      } else {# traits is a trait, insert in a list
        traits <- list(traits)
      }

      traitsNames <- vapply(traits, function(t){t$name}, "name")
      if (!identical(traitsNames, unique(traitsNames))) {
        stop('All traits names should be differents. Please check the the "name" field of your traits')
      }
      names(traits) <- traitsNames



      if (!is.numeric(plotCost) || length(plotCost) != 1 || plotCost < 0) {
        stop("plotCost should be a single positive numeric value")
      }

      if (!is.numeric(mu)) {
        stop('"mu" should be a single numeric value or a numeric vector of length equal to the number of traits')
      } else if (length(mu) == 1) {
        mu <- rep(mu, length(traits))
      } else if (length(mu) != length(traits)) {
        stop('"mu" should be a singler numeric value or a numeric vector of length equal to the number of traits')
      }
      if (is.null(names(mu))) {
        names(mu) <- traitsNames
      } else if (!identical(sort(names(mu)), sort(traitsNames))) {
        stop('names(mu) should be compatible with the traits names')
      } else {
        mu <- mu[traitsNames]
      }


      if (is.null(ve) && is.null(he) && is.null(pop)) {
        stop('you should either provide "ve" or "he" and "pop"')
      }

      if (is.numeric(ve)) {

        if (!is.null(he) || !is.null(pop) ) {
          stop('If "ve" is specified, "he" and "pop" should not be specify.')
        }

        if (length(ve) == 1) {
          ve <- rep(ve, length(traits))
        } else if (length(ve) != length(traits)) {
          stop('length(ve) should be equal to 1 or to the number of traits.')
        }

        if (is.null(names(ve))) {
          names(ve) <- traitsNames
        } else if (!identical(sort(names(ve)), sort(traitsNames))) {
          stop('names(ve) should be compatible with the traits names')
        } else {
          ve <- ve[traitsNames]
        }

      } else if (is.null(ve)) {
        if (is.null(he) || is.null(pop)) {
          stop('If you do not provide "ve" please specify both "he" and "pop"')
        }

        if (any(he > 1) || any(he < 0)) {
          stop('"he" should be between 0 and 1')
        }

        if (length(he) == 1) {
          he <- rep(he, length(traits))
        } else if (length(he) != length(traits)) {
          stop('length(he) should be equal to 1 or to the number of traits.')
        }

        if (is.null(names(he))) {
          names(he) <- traitsNames
        } else if (!identical(sort(names(he)), sort(traitsNames))) {
          stop('names(he) should be compatible with the traits names')
        } else {
          he <- he[traitsNames]
        }

        if (class(pop)[1] != "population") {
          stop('pop should be an object of class "population"')
        }
      } else {
        stop('Please provide "ve" or "he" and "pop"')
      }


      self$name <- name

      self$traits <- traits
      private$traitsNames <- traitsNames

      self$plotCost <- plotCost
      self$mu <- mu


      # Calculation of ve
      if (is.numeric(ve)) {
        self$ve <- ve
      } else {
        ve <- vapply(traits, function(trait) {
          g <- trait$gv(pop)
          ve <- var(g) * ((1 / he[trait$name]) - 1)
          ve
        }, 1)
        names(ve) <- traitsNames
        self$ve <- ve
      }


    },
    #' @description
    #' Display informations about the object
    print = function() {
      cat(paste0("Phenotyper: ", self$name, "\n",
                 "Traits: ", paste(private$traitsNames, collapse = ", "), "\n",
                 "\u03bc: ", paste(self$mu, collapse = ", "), "\n",
                 "\u03c32: ", paste(self$ve, collapse = ", "), "\n",
                 "Phenotyping cost (per plot): ", self$plotCost, "\n",
                 "timeEffect: ", self$timeEffect))
    },
    #' @description Phenotype a given population
    #' @param pop [population class] population
    #'     (see: \link[breedSimulatR]{population})
    #' @param rep [numeric] number of replication for each individuals
    #' @param offset [numeric] offset added to the phenotypic calculation
    #' results
    #' @details
    #' phenotypic values for individual i repetition j is calculated as follow:
    #' \deqn{y_{ij} = \mu + g_i + e_ij + offset}
    #' \deqn{e_ij \sim N(0, \sigma_e^2)}
    #' where g_i is the genetic value of the individual i
    #' @examples
    #'
    #' pheno1 <- phenoLab1$trial(example_pop, rep = 4, offset = c(3, 5))
    #' pheno1$cost
    #' summary(pheno1$data)
    #'
    #' pheno2 <- phenoLab2$trial(example_pop, rep = 3, offset = c(-5, 0))
    #' pheno2$cost
    #' summary(pheno2$data)
    #'
    #' pheno3 <- phenoLab2$trial(example_pop,
    #'                           rep = round(runif(example_pop$nInd, 1, 3)),
    #'                           offset = c(-5, 0))
    #' pheno3$cost
    #' summary(pheno3$data)
    trial = function(pop, rep = 1, offset = 0) {

      if (class(pop)[1] != "population") {
        stop('pop should be an object of class "population"')
      }

      if (!is.numeric(rep)) {
        stop('rep should be a numeric vector')
      }
      if (length(rep) != 1 && length(rep) != pop$nInd) {
        stop('length(rep) should be equal to 1 or to the number of
        individuals in the population')
      }
      if (length(rep) == 1) {
        rep <- base::rep(rep, pop$nInd)
      }

      if (any(rep < 0)) {
        stop('all values of rep should be higher or equal than 0')
      }

      if (!is.numeric(offset)) {
        stop('"offset" should be a numeric value')
      }
      if (length(offset) == 1) {
        offset <- base::rep(offset, length(self$traits))
      } else if (length(offset) != length(self$traits)) {
        stop('"length(offset)" should be equal to the number of traits')
      }
      if (is.null(names(offset))) {
        names(offset) <- private$traitsNames
      } else if (!identical(sort(names(offset)), sort(private$traitsNames))) {
        stop('names(offset) should be compatible with the traits names')
      } else {
        offset <- offset[private$traitsNames]
      }

      cost <- self$plotCost * sum(rep)
      pheno <- matrix(sapply(self$traits, function(trait) {
        mu <- self$mu[trait$name]
        sigma <- sqrt(self$ve[trait$name])

        gv <- base::rep(trait$gv(pop), rep)
        e <- rnorm(sum(rep), sd = sigma)

        matrix(mu + gv + e + offset[trait$name], ncol = 1)
      }),
      nrow = sum(rep))
      colnames(pheno) <- names(self$traits)
      dta <- data.frame(ind = base::rep(row.names(pop$genoMat), rep),
                        pheno,
                        rep = unlist(lapply(rep, function(x){
                          if (x != 0) {
                            out <- seq(x)
                          } else {
                            out <- c()
                          }
                          out
                        })),
                        phenotyper = self$name)

      list(
        data = dta,
        cost = cost
      )



    },
    #' @description Calculate the heritability of a given population
    #' @param pop [population] the population
    #' @return named numeric vector of the calculated heritability of the
    #' population for each trait
    #' @examples
    #' phenoLab1$he(example_pop)
    #' phenoLab2$he(example_pop)
    he = function(pop) {
      sapply(self$traits, function(t){
        g <- t$gv(pop)
        he <- var(g) / (var(g) + self$ve[t$name])
        he
      })
    }
  ),
  private = list(
    ## Private Fields ####
    # @field traitsNames [character] names of the traits
    traitsNames = NULL
  )
)
