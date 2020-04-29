# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# Definition of Population class



#' R6 class representing a population
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
    #' rawHaplo1 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
    #'                     nrow = 2)
    #' colnames(rawHaplo1) <- sprintf(fmt = paste0("SNP%0", 2,"i"),
    #'                               1:(3 + 4 + 5))
    #' haplo1 <- haplotype$new(SNPinfo = SNPs,
    #'                        haplo = rawHaplo1)
    #' rawHaplo2 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
    #'                     nrow = 2)
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
    initialize = function(name = NULL,
                          inds = list(),
                          verbose = TRUE){
      # checks
      if (is.null(name)) {
        name <- "Unspecified"
      } else name <- as.character(name)

      if (class(inds)[1] == "individual") {
        inds <- list(inds)
      } else if (class(inds) != "list") {
        stop(paste("inds must be an individual object or a list of",
                   "individuals objects"))
      }

      if (verbose) {
        cat("Create population: Add individuals...\n")
        i <- 1
        tot <- length(inds)
      }
      self$name <- name
      for (ind in inds) {
        if (verbose) {
          cat(paste0("\r", round(i/tot*100), "%"))
          i <- i + 1
        }
        private$addInd(ind)
      }
      if (verbose) cat(paste("\nA new population created:", self$name, "!\n"))
    },
    #' @description Add individuals to the population
    #' @param inds [individual class or list] list of individuals of the
    #'     population (see:\link[breedSimulatR]{individual})
    #' @examples
    #' # create new individual
    #'
    #' rawHaplo3 <- matrix(sample(c(0, 1), (3 + 4 + 5) * 2, replace = TRUE),
    #'                     nrow = 2)
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
    #' @param indsNames [character] character vetcor of the individuals' names
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
    },
    #' @field genoMat [matrix] matrix of all the genotypes of the population
    #'   encoded in allele doses. (individuals in row and markers in column)
    genoMat = function(){
      if (length(self$inds) > 0) {
        t(vapply(self$inds, function(ind){
          ind$haplo$allelDose
        }, vector(mode = "numeric",
                  length = length(self$inds[[1]]$haplo$allelDose)))
        )
      } else {
        NULL
      }

    },

    #' @field af [named vector] allele frequency
    af = function(){
      ploidy <- self$specie$ploidy
      stopifnot(ploidy %in% c(1,2))

      geno <- self$genoMat
      freq <- colSums(geno) / (self$nInd * ploidy)
      freq
    },

    #' @field maf [named vector] minor allele frequency
    maf = function(){
      ploidy <- self$specie$ploidy
      stopifnot(ploidy %in% c(1,2))
      freq <- 0.5 - abs(self$af - 0.5)
      freq
    }

  ),
  private = list(
    # @description Add new individual to the population
    # @param ind [individual class] individual
    #   (see:\link[breedSimulatR]{individual})
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
        stop(paste("Individual of a different species than the population's",
                   "one.\nPlease add", self$specie$name, "individuals."))
      }

      # check name
      if (ind$name %in% names(self$inds)) {
        stop(paste("Individual with the same name already exists",
                   "in the population:", ind$name))
      }

      # add new individual
      self$inds[[ind$name]] <- ind
      NULL
    }
  )
)





#' Create population object from genotype data.frame
#'
#' @param geno [data.frame] genotype of the individuals encoded in allele dose.
#'   All individuals should be homozygotes. (value 0 or 2)
#' @param SNPinfo [SNPinfo object] information about the individuals haplotypes'
#'   SNPs (see:\link[breedSimulatR]{SNPinfo})
#' @param indNames NULL or character string vector specifying the individuals
#'   names. If NULL, \code{rownames(geno)} will be used.
#' @param popName [character string] population's name.
#' @param verbose [boolean] display information
#' @return population object (see:\link[breedSimulatR]{population})
#' @export
#'
#' @examples
#' mySpec <- specie$new(nChr = 10,
#'                      lchr = 10^6,
#'                      specName = "Geneticae Exempli",
#'                      ploidy = 2,
#'                      recombRate = 3/10^6)
#' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
#'                     specie = mySpec)
#'
#' print(exampleData$genotypes[1:5, 1:10])
#' example_pop <- createPop(geno = exampleData$genotypes,
#'                          SNPinfo = SNPs,
#'                          popName = "Example population")
createPop <- function(geno,
                      SNPinfo,
                      indNames = NULL,
                      popName = NULL,
                      verbose = TRUE) {
  if (verbose) {
    cat("Create population: Initialisation...\n")
  }
  # check parameters:
  if (any(!colnames(geno) %in% SNPinfo$SNPcoord$SNPid)) {
    stop('Some markers of "geno" are not in "SNPinfo"')
  }

  if (any(!SNPinfo$SNPcoord$SNPid %in% colnames(geno))) {
    warning('Some markers of "SNPinfo" are not in "geno"')
  }

  if (is.null(indNames)) {
    if (is.null(rownames(geno))) {
      stop('"rownames(geno)" is NULL, please specify "indNames"')
    } else indNames <- rownames(geno)
  } else if (length(indNames) == 1) {
    indNames <- sprintf(
      fmt = paste0(indNames,
                   "%0", floor(log10(nrow(geno))) + 1, "i"),
      seq(nrow(geno)))
  } else if (length(indNames) != nrow(geno)) {
    stop(paste0("length(indNames) = ", length(indNames),
                '\n"length(indNames)" must be equal to "1" or "nrow(geno)"'))
  }

  if (!all.equal(unique(indNames), indNames)) {
    stop('All values of "indNames" must be different')
  }

  if (!all(unique(unlist(geno)) %in% c(0,2))) {
    stop(paste0('Some values of "geno" are different from 0 or 2.\n',
                'All individulas must be homozygotes and genotypes must be ',
                'encoded as allele doses'))
  }

  listInds <- list()
  geno <- as.matrix(geno)

  if (verbose) {
    cat("Create population: Create individuals...\n")
    prog <- 0
    nInd <- nrow(geno)
  }

  for (i in seq(nrow(geno))) {
    if (verbose) {
      prog <- i/nInd
      cat(paste0("\r", round(prog*100), "%"))
    }

    haplo <- geno[i,]
    listInds[[i]] <- individual$new(
      name = indNames[i],
      specie = SNPinfo$specie,
      parent1 = NA,
      parent2 = NA,
      haplo = haplotype$new(SNPinfo, rbind(haplo, haplo) / 2),
      verbose = FALSE
    )
  }

  if (verbose) {
    cat("\nCreate population: Create population object...\n")
    prog <- 0
  }
  population$new(name = popName, inds = listInds, verbose = verbose)

}
