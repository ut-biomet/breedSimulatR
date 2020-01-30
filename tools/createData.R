# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# R script to generate package's data

library(dplyr)
library(breedSimulatR)

# load genotype data:
geno <- read.table("tools/data/genotypes.txt.gz", header = T, sep = "\t")

# load snp coordinates:
snpCoord <- read.table("tools/data/snp_coords.txt.gz", header = T, sep = "\t")
snpCoord$SNPid <- row.names(snpCoord)
snpCoord$chr <- gsub("c", "C", snpCoord$chr)
snpCoord$chr <- gsub("Chr(?=\\d$)", "Chr0", snpCoord$chr, perl = TRUE)


# check homogeneity between geno and snp cood
if (any(!colnames(geno) %in% snpCoord$SNPid)) {
  stop('Some markers of "geno" are not in "snpCoord".')
}

if (any(!snpCoord$SNPid %in% colnames(geno))) {
  warning('Some markers of "snpCoord" are not in "geno".')
}


# create specie
snpCoord %>%
  group_by(chr) %>%
  summarise(len = max(pos))

example_Specie <- specie$new(specName = "Statisticae exempli",
                             nChr = 10,
                             lchr = 1000000,
                             ploidy = 2,
                             recombRate = 3/1000000)

example_SNPs <- SNPinfo$new(SNPcoord = snpCoord,
                            specie = example_Specie)

example_pop <- createPop(geno = geno,
                         SNPinfo = example_SNPs,
                         popName = "Example population")

example_genotypes <- geno


all.equal(as.matrix(geno[,sort(colnames(geno))]), example_pop$genoMat)
colnames(geno)[1:5]
colnames(example_pop$genoMat)[1:5]
length(example_pop$inds$Coll0001$haplo$allelDose)


# usethis::use_data(example_Specie, overwrite = TRUE)
# usethis::use_data(example_SNPs, overwrite = TRUE)
# usethis::use_data(example_pop, overwrite = TRUE)
# usethis::use_data(example_genotypes, overwrite = TRUE)
