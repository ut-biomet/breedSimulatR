% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population.R
\name{createPop}
\alias{createPop}
\title{Create population object from genotype data.frame}
\usage{
createPop(geno, SNPinfo, indNames = NULL, popName = NULL, verbose = TRUE)
}
\arguments{
\item{geno}{[data.frame] genotype of the individuals encoded in allele dose.
All individuals should be homozygotes. (value 0 or 2)}

\item{SNPinfo}{[SNPinfo object] information about the individuals haplotypes'
SNPs (see:\link[breedSimulatR]{SNPinfo})}

\item{indNames}{NULL or character string vector specifying the individuals
names. If NULL, \code{rownames(geno)} will be used.}

\item{popName}{[character string] population's name.}

\item{verbose}{[boolean] display information}
}
\value{
population object (see:\link[breedSimulatR]{population})
}
\description{
Create population object from genotype data.frame
}
\examples{
mySpec <- specie$new(nChr = 10,
                     lchr = 10^6,
                     lchrCm = 100,
                     specName = "Geneticae Exempli")
SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
                    specie = mySpec)

print(exampleData$genotypes[1:5, 1:10])
example_pop <- createPop(geno = exampleData$genotypes,
                         SNPinfo = SNPs,
                         popName = "Example population")
}
