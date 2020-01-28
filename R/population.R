# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of Population class



#' An R6 class representing a population
#'
#' @description
#' population object store specific information about one individual
#'
#'
#' @export
#' @import R6
population <- R6::R6Class(
  "population",
  lock_objects = FALSE,
  public = list(
    #' @field name [string] Name of the population
    name = NULL,
    #' @field specie [specie class] Specie of the SNPs
    #'   (see:\link[breedSimulatR]{specie})
    specie = NULL,
    #' @field inds [list] list of population's individuals
    inds = list(),
    #' @description Create a new population object.
    #' @param name [string] name of the population
    #' @param inds [individual class or list] list of individuals of the
    #'     population (see:\link[breedSimulatR]{individual})
    #' @param verbose [boolean] display information
    #' @return A new `population` object.
    #' @examples
    #' # create specie:
    #' mySpec <- specie$new(nChr = 3,
    #'                      lchr = c(100, 150, 200),
    #'                      recombRate = 0.006,
    #'                      verbose = FALSE)
    #'
    #' # simulate SNP:
    #' SNPcoord <- data.frame(chr = c(rep("Chr1", 3),
    #'                                rep("Chr2", 4),
    #'                                rep("Chr3", 5)),
    #'                        pos = c(sample(100, 3),
    #'                                sample(150, 4),
    #'                                sample(200, 5)),
    #'                        SNPid = sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                                        1:(3 + 4 + 5)))
    #'
    #' # create SNPinfo object:
    #' SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = mySpec)
    #'
    #'
    #' # simulate haplotype:
    #' rawHaplo1 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE), nrow = 2)
    #' colnames(rawHaplo1) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                               1:(3 + 4 + 5))
    #' haplo1 <- haplotype$new(SNPinfo = SNPs,
    #'                        haplo = rawHaplo1)
    #' rawHaplo2 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE), nrow = 2)
    #' colnames(rawHaplo2) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                               1:(3 + 4 + 5))
    #' haplo2 <- haplotype$new(SNPinfo = SNPs,
    #'                        haplo = rawHaplo2)
    #'
    #'
    #' # create individuals:
    #' myInd1 <-  individual$new(name = "Ind 1",
    #'                          specie = mySpec,
    #'                          parent1 = "OkaaSan",
    #'                          parent2 = "OtouSan",
    #'                          haplo = haplo1,
    #'                          verbose = FALSE)
    #' myInd2 <-  individual$new(name = "Ind 2",
    #'                          specie = mySpec,
    #'                          parent1 = "OkaaSan",
    #'                          parent2 = "OtouSan",
    #'                          haplo = haplo2,
    #'                          verbose = FALSE)
    #' myPop <- population$new(name = "My Population 1",
    #'                         inds = list(myInd1, myInd2),
    #'                         verbose = FALSE)
    initialize = function(name = "New Pop",
                          inds = list(),
                          verbose = T){
      # checks
      name <- as.character(name)

      if (class(inds)[1] == "individual") {
        inds <- list(inds)
      } else if (class(inds) != "list") {
        stop("inds must be an individual object or a list of individuals objects")
      }

      self$name = name
      for (ind in inds) {
        private$addInd(ind)
      }
      if (verbose) cat(paste("A new population created:", self$name, "!"))
    },
    #' @description Add individuals to the population
    #' @param inds [individual class or list] list of individuals of the
    #'     population (see:\link[breedSimulatR]{individual})
    #' @examples
    #' # create new individual
    #'
    #' rawHaplo3 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE), nrow = 2)
    #' colnames(rawHaplo3) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                               1:(3 + 4 + 5))
    #' haplo3 <- haplotype$new(SNPinfo = SNPs,
    #'                        haplo = rawHaplo3)
    #' myInd3 <-  individual$new(name = "Ind 3",
    #'                          specie = mySpec,
    #'                          parent1 = "OkaaSan",
    #'                          parent2 = "OtouSan",
    #'                          haplo = haplo1,
    #'                          verbose = FALSE)
    #'
    #' # add individual
    #' print(myPop)
    #' myPop$addInds(myInd3)
    #' print(myPop)
    addInds = function(inds){
      # checks
      if (class(inds)[1] == "individual") {
        inds <- list(inds)
      }
      lapply(inds, private$addInd)

      invisible(NULL)

    },
    #' @description Remove individuals from the population
    #' @param indsNames [character] character vercor of the individuals' names
    #' @examples
    #' print(myPop)
    #' myPop$remInds("Ind 2")
    #' print(myPop)
    remInds = function(indsNames){

      # check
      if (class(indsNames) != "character") {
        stop("Please provide individuals' names as a character vector.")
      }

      if (any(! indsNames %in% names(self$inds))) {
        id <- which(! indsNames %in% names(self$inds))
        warning(paste("Some individuals to remove are not in the population:",
                      paste(indsNames[id], collapse = " ; ")))
      }

      self$inds[indsNames] <- NULL
      invisible(NULL)
    },
    #' @description
    #' Display informations about the object
    print = function() {
      cat(paste0(
        "Population: ", self$name, "\n",
        "Species: ", self$specie$specName, "\n",
        "Number of individuals: ", self$nInd
      ))
    }


  ),
  active = list(
    #' @field nInd [numeric] number of individual in the population
    nInd = function(){
      length(self$inds)
    }

  ),
  private = list(
    # @description Add new individual to the population
    # @param ind [individual class] individual (see:\link[breedSimulatR]{individual})
    # @return NULL
    addInd = function(ind) {

      # checks class
      if (class(ind)[1] != "individual") {
        stop('variable "ind" must be of class "individual".')
      }

      # check species
      if (is.null(self$specie)) {
        self$specie <- ind$specie
      } else if (!isTRUE(all.equal(self$specie, ind$specie))) {
        stop(paste("Individual of a different species than the population's one.\n",
                   "Please add", self$specie$name, "individuals."))
      }

      # check name
      if (ind$name %in% names(self$inds)) {
        stop(paste("Individual with the same name already exists in the population:", ind$name))
      }

      # add new individual
      self$inds[[ind$name]] <- ind
      NULL
    }
  )
)
