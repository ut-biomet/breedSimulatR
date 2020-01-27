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
    #' @description Create a new population object.
    #' @param name [string] name of the population
    #' @param inds [individual class or list] list of individuals of the
    #'     population (see:\link[breedSimulatR]{individual})
    #' @param verbose [boolean] display information
    #' @return A new `population` object.
    #' @examples
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
    addInds = function(inds, verbose = T){
      # checks
      if (class(inds)[1] == "individual") {
        inds <- list(inds)
      }

      lapply(inds, private$addInd)
      NULL

    },
    #' @description Remove individuals from the population
    #' @param indsNames [character] character vercor of the individuals' names
    #' @examples
    remInds = function(indsNames, verbose = T){

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

    },
    #' @description
    #' Display information about the object
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
    #' @description Add new individual to the population
    #' @param ind [individual class] individual (see:\link[breedSimulatR]{individual})
    #' @return NULL
    #' @examples
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
