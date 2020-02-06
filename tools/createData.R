# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# R script to generate package's data

library(dplyr)
library(breedSimulatR)

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
                             ploidy = 2,
                             recombRate = 3/1e6)

# create SNPinfo
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
