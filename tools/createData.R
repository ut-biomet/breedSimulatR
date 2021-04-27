# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# R script to generate package's data

library(dplyr)
# library(breedSimulatR)
devtools::load_all()

# load genotype data:
example_genotypes <- read.table("tools/data/genotypes.txt.gz", header = T, sep = "\t")

# load snp coordinates:
example_snpCoord <- read.table("tools/data/snp_coords.txt.gz", header = T, sep = "\t")

# load snp effects:
snpEffects <- read.table("tools/data/snpEffects.txt", header = T, sep = "\t")
example_snpEffects <- snpEffects$x
names(example_snpEffects) <- rownames(snpEffects)

# check homogeneity between geno and snp cood
if (any(!colnames(example_genotypes) %in% example_snpCoord$SNPid)) {
  stop('Some markers of "example_genotypes" are not in "example_snpCoord".')
}

if (any(!example_snpCoord$SNPid %in% colnames(example_genotypes))) {
  warning('Some markers of "example_snpCoord" are not in "example_genotypes".')
}


# create specie
example_specie <- specie$new(specName = "Statisticae exempli",
                             nChr = 10,
                             lchr = 1e6,
                             lchrCm = 100)

# create SNPinfo

# simulate linkage map position
b1 <- c(6, 13, 10, 7, 13, 10, 13, 9, 6, 9)
b2 <- c(6, 8, 12, 12, 11, 7, 9, 7, 7, 13)
names(b1) <- names(b2) <- example_specie$chrNames
example_snpCoord <- do.call(rbind,
  lapply(example_specie$chrNames,
         function(chr){
           SNPcoord <- example_snpCoord[example_snpCoord$chr == chr,]
           SNPcoord$linkMapPos <- breedSimulatR:::.simulLinkMapPos(
             SNPcoord$physPos,
             example_specie$lchr[chr],
             example_specie$lchrCm[chr],
             b1[chr],
             b2[chr])
           SNPcoord[c("chr", "physPos", "linkMapPos", "SNPid")]
         }))


example_SNPs <- SNPinfo$new(SNPcoord = example_snpCoord,
                            specie = example_specie)

# create pop
example_pop <- createPop(geno = example_genotypes,
                         SNPinfo = example_SNPs,
                         popName = "Example population")


# checks
stopifnot(all.equal(as.matrix(example_genotypes[,sort(colnames(example_genotypes))]), example_pop$genoMat))
stopifnot(all.equal(names(example_snpEffects), colnames(example_pop$genoMat)))



exampleData <- list(
  snpCoord = example_snpCoord,
  genotypes = example_genotypes,
  snpEffects = example_snpEffects
)


file.remove(list.files("data/", full.names = TRUE))
usethis::use_data(exampleData, overwrite = TRUE)
