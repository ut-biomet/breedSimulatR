# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# File defining usefull function for test scripts



##### Initialisation functions ####
create_spec <- function(nChr = round(runif(1, 1, 10)),
                        lchr = round(pmax(rnorm(nChr, 450, 50), 200)),# > 200
                        lchrCm = 100,
                        name = "Undefinded") {

  specie$new(nChr = nChr,
             lchr = lchr,
             lchrCm = lchrCm,
             specName = name,
             verbose = F)
}

create_SNP <- function(spec = create_spec(), nMarker = NULL,
                       lmap = TRUE, b1 = runif(1,5,15), b2 = runif(1,5,15)){

  if (is.null(nMarker)) {
    nMarker <- round(spec$lchr/10)
  } else {
    stopifnot(all(nMarker < spec$lchr))
    nMarker <- rep(nMarker, spec$nChr)
  }

  # generate arbitrary marker position
  physPos <- c()
  for (i in seq(spec$nChr)) {
    physPos <- c(physPos, sample(spec$lchr[i], nMarker[i]))
  }



  SNPcoord <- do.call(rbind,lapply(seq(spec$nChr),
                                   function(chr){
                                     physPos <- sample(spec$lchr[chr], nMarker[chr],
                                                       replace = FALSE)
                                     if (lmap) {
                                       linkMapPos <- .simulLinkMapPos(physPos,
                                                                      spec$lchr[chr],
                                                                      spec$lchrCm[chr],
                                                                      b1 = b1,
                                                                      b2 = b2)
                                     } else {
                                       linkMapPos <- NA
                                     }

                                     data.frame(physPos,linkMapPos)
                                   }))
  SNPcoord$chr <- rep(spec$chrNames, times = nMarker)
  # generate arbitrary SNPid
  SNPcoord$SNPid <- .charSeq("SNP", sample(sum(nMarker)*50, sum(nMarker)))

  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = spec)
  SNPs
}

create_haplo <- function(SNPs, af = NULL){

  if (is.null(af)) {
    af <- 0.5
  } else {
    stopifnot(length(af) == SNPs$nSNP())
  }

  rawHaplo <- matrix(rbinom(SNPs$nSNP() * SNPs$specie$ploidy, 1,
                            af),
                     nrow = SNPs$specie$ploidy)
  colnames(rawHaplo) <- SNPs$SNPcoord$SNPid
  haplo <- haplotype$new(SNPinfo = SNPs,
                         haplo = rawHaplo)
}


create_inds <- function(haploList){
  if (class(haploList)[1]  == "Haplotype") {
    haploList <- list(haploList)
  }

  spec <- haploList[[1]]$SNPinfo$specie
  inds <- mapply(function(haplo, id){
    myInd <- individual$new(name = paste("Ind", id),
                            specie = spec,
                            parent1 = paste("OkaaSan", id),
                            parent2 = paste("OtouSan", id),
                            haplo = haplo,
                            verbose = F)
  },
  haploList, c(seq_along(haploList)))

  if (length(inds) == 1) {
    return(inds[[1]])
  }
  inds
}


create_trait <- function(SNPs,
                         name = "trait",
                         n = 100,
                         qtn = NULL,
                         qtnEff = NULL) {
  if (is.null(qtn)) {
    qtn <-  sample(SNPs$SNPcoord$SNPid, n)
  }

  if (is.null(qtnEff)) {
    qtnEff <-  rnorm(length(qtn), sd = 0.5)
  }

  myTrait <- trait$new(
    name = name,
    qtn = qtn,
    qtnEff = qtnEff
  )
  myTrait
}

.charSeq <- breedSimulatR:::.charSeq
.simulLinkMapPos <- breedSimulatR:::.simulLinkMapPos
