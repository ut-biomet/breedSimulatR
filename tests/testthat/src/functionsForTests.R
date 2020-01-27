# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 / 2020 The University of Tokyo
#
# Description:
# File defining usefull function for test scripts



##### Initialisation functions ####
create_spec <- function(nChr = round(runif(1, 1, 10)),
                        lchr = round(pmax(rnorm(nChr, 450, 50), 200)),# at least 200 markers
                        ploidy = 2,
                        recombRate = 3/sum(lchr),
                        name = "Undefinded"){
  specie$new(nChr = nChr,
             lchr = lchr,
             ploidy = ploidy,
             recombRate = recombRate,
             specName = name,
             verbose = F)
}

create_SNP <- function(spec){

  nMarker <- round(spec$lchr/10)

  # generate arbitrary marker position
  pos <- c()
  for (i in seq(spec$nChr)) {
    pos <- c(pos, sample(spec$lchr[i], nMarker[i]))
  }

  # generate arbitrary SNPid
  SNPid <- sprintf(fmt = paste0("SNP%0", ceiling(log10(sum(nMarker)*50)),"i"),
                   sample(sum(nMarker)*50, sum(nMarker)))
  SNPcoord <- data.frame(chr = rep(spec$chrNames, times = nMarker),
                         pos = pos,
                         SNPid = SNPid)
  # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = SNPcoord, specie = spec)
  SNPs
}

create_haplo <- function(SNPs){
  rawHaplo <- matrix(sample(c(0, 1), SNPs$nSNP() * SNPs$specie$ploidy, replace = T),
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
  haploList, c(1:length(haploList)))

  if (length(inds) == 1) {
    return(inds[[1]])
  }
  inds
}