# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Definition of functions reading external data file and creating
# breedSimulatR objects

#' Create `population` and `SNPinfo` objects from VCF file
#'
#' @param file path of a `.vcf` or `.vcf.gz` file.
#' @param specie [specie class] Specie of the population
#'   (see:\link[breedSimulatR]{specie}), if this parameter is not provided the function will create one based on the information in the VCF file. See details section for more information.
#' @param verbose If TRUE, display information on the function progress
#'
#' @details
#' If `specie` is not provided, the function will create one based on the information in the VCF file. The number of chromosome will be set to the number of differents chromosome in the VCF file, the names of the chromosomes will be set the the unique values of the column `CHROM` of the VCF file, the physical length of each chromosome will equal to the highest physical position of the chromosome's markers, the length of chromosomes in centimorgans will be set to 1cm per Mega Base.
#' @return list of 3 elements `pop` of class `population`, `snps` of class `snpInfo` and `specie` of class specie
#' @export
#'
#' @importFrom vcfR read.vcfR
#' @examples
#' \dontrun{
#' # create vcf.gz file
#' mySpec <- specie$new(nChr = 10,
#'                      lchr = 10^6,
#'                      lchrCm = 100,
#'                      specName = "Geneticae Exempli")
#' SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
#'                     specie = mySpec)
#' example_pop <- createPop(geno = exampleData$genotypes,
#'                          SNPinfo = SNPs,
#'                          popName = "Example population")
#' vcfFile <- tempfile(fileext = "vcf.gz") 
#' example_pop$writeVcf(vcfFile)
#'
#' # read vcf file
#' x <- readVCF(file = vcfFile , verbose = TRUE)
#' print(x$pop)
#' print(x$snps)
#' print(x$specie)
#' }
readVCF <- function(file, specie = NULL, verbose = TRUE) {

  # checks
  if (!is.null(specie) && class(specie)[1] != "Specie") {
    stop('"class(specie)" must be "Specie"')
  }


  # read vcf file
  if (verbose) {
    cat("Read VCF file using `vcfR` package\n")
  }
  # browser()
  vcf <- vcfR::read.vcfR(file = file, verbose = verbose)
  fixVcf <- as.data.frame(vcfR::getFIX(vcf))
  fixVcf[,"POS"] <- as.numeric(fixVcf[,"POS"])

  # create specie
  if (is.null(specie)) {
    if (verbose) {
      cat("`specie` is not provided, creation of a new specie based on the information in the vcf file\n")
    }
    chrNames <- unique(fixVcf[,"CHROM"])
    nChr <- length(chrNames)

    # get chromosome length
    tabLchr <- stats::aggregate(POS ~ CHROM, data = fixVcf, max)
    lchr <- c()
    lchr[tabLchr[,"CHROM"]] <- tabLchr[,"POS"]
    lchr <- lchr[chrNames]

    specie <- breedSimulatR::specie$new(nChr = nChr,
                         chrNames = chrNames,
                         lchr = lchr,
                         lchrCm = lchr/10^6, # 1 cMper Mb.
                         verbose = verbose)
  }



  # extract SNP informations
  SNPcoord <- fixVcf[,c("CHROM", "POS", "ID")]
  colnames(SNPcoord) <- c("chr", "physPos", "SNPid")
  SNPcoord$linkMapPos <- NA
  SNPcoord <- SNPcoord[,c("chr", "physPos", "linkMapPos", "SNPid")]

  snps <- SNPinfo$new(SNPcoord = SNPcoord, specie = specie)



  # extract haplotypes
  # haps <- vcfR::extract.haps(vcf, return.alleles = FALSE)
  if (verbose) {
    cat("Exctract genetic information...\n")
  }
  if (all(vcf@gt[,"FORMAT"] == "GT")) {
    # quick return if FORMAT == GT
    gt <- vcf@gt
    gt <- gt[, colnames(gt) != "FORMAT"]
    row.names(gt) <- vcf@fix[,"ID"]
  } else {
    gt <- vcfR::extract.gt(vcf, element = 'GT')
  }

  # stop if not phased
  if (!all(grepl("|", gt, fixed = TRUE))) {
    stop("VCF file should be phased for all variant and all individuals, (`|` separator for GT filed).")
  }


  rm(vcf) # free memory
  if (verbose) {
    cat("Exctract haplotypes... \n")
  }
  indNames <- colnames(gt)
  markersNames <- rownames(gt)
  gt <- paste(gt, collapse = "|") # make strsplit faster
  gt <- strsplit(gt, split = "|", fixed = TRUE) # faster when using fixed=T
  gt <- unlist(gt)
  gt <- as.integer(gt)
  # values of gt are mixed so we need to reorder:
  gt <- c(gt[seq(1, length(gt), 2)], gt[seq(2, length(gt), 2)])
  gt <- matrix(gt,
               nrow = length(markersNames),
               byrow = FALSE)
  rownames(gt) <- markersNames
  colnames(gt) <- paste(rep(indNames, 2),
                         rep(1:2, each=length(indNames)),
                         sep = "_")

  # create list of individuals
  if (verbose) {
    cat("Create individuals...\n")
  }
  listInds <- vector(mode = "list", length = length(indNames))
  names(listInds) <- indNames
  i <- 1
  for (indName in indNames) {
    if (verbose) {
      prog <- i/length(indNames)
      cat(paste0("\r", round(prog*100), "%"))
      i <- i+1
    }
    # browser()
    haplo <- t(gt[, paste(indName, c(1,2), sep = "_")])
    haplo <- haplotype$new(SNPinfo = snps,
                           haplo = haplo)
    listInds[[indName]] <- individual$new(name = indName,
                                          specie = specie,
                                          haplo = haplo,
                                          verbose = FALSE)
  }
  if (verbose) {cat("\n")}

  # create population
  pop <- population$new(name = "", inds = listInds, verbose = verbose)


  # return results
  out <- list(pop = pop,
              snps = snps,
              specie = specie)
  return(out)
}
