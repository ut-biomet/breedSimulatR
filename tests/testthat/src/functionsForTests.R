# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# File defining usefull function for test scripts



##### Initialisation functions ####
create_spec <- function(nChr = round(runif(1, 1, 10)),
                        lchr = round(pmax(rnorm(nChr, 450, 50), 200)),# > 200
                        ploidy = 2,
                        recombRate = 3 / sum(lchr),
                        name = "Undefinded") {

  specie$new(nChr = nChr,
             lchr = lchr,
             ploidy = ploidy,
             recombRate = recombRate,
             specName = name,
             verbose = F)
}

create_SNP <- function(spec, nMarker = NULL){

  if (is.null(nMarker)) {
    nMarker <- round(spec$lchr/10)
  } else {
    stopifnot(all(nMarker < spec$lchr))
    nMarker <- rep(nMarker, spec$nChr)
  }

  # generate arbitrary marker position
  pos <- c()
  for (i in seq(spec$nChr)) {
    pos <- c(pos, sample(spec$lchr[i], nMarker[i]))
  }

  # generate arbitrary SNPid
  SNPid <- .charSeq("SNP", sample(sum(nMarker)*50, sum(nMarker)))
  SNPcoord <- data.frame(chr = rep(spec$chrNames, times = nMarker),
                         pos = pos,
                         SNPid = SNPid)
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

.charSeq <- breedSimulatR:::.charSeq
