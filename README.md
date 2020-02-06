
<!-- README.md is generated from README.Rmd. Please edit that file -->

# breedSimulatR

<!-- badges DEV : start -->

[![pipeline
status](https://gitlab.com/juliendiot42/breedSimulatR/badges/dev/pipeline.svg)](https://gitlab.com/juliendiot42/breedSimulatR/commits/dev)
[![codecov](https://codecov.io/gh/ut-biomet/breedSimulatR/branch/dev/graph/badge.svg?token=4Uxp1ySLIn)](https://codecov.io/gh/ut-biomet/breedSimulatR)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN
status](https://www.r-pkg.org/badges/version/breedSimulatR)](https://CRAN.R-project.org/package=breedSimulatR)
<!-- badges: end -->

R package providing classes and functions to simulate breeding schemes.

## Installation

:octocat: You can install `breedSimulatR` from
[GitHub](https://github.com/) with:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("ut-biomet/breedSimulatR")
```

You can check the installation with these lines:

``` r
library(breedSimulatR)
help(package = breedSimulatR)
```

## Example

This is a basic example which shows how to use the package.

The package contains some example data that we will use for this
example. These data are stored in the varaible `exampleData`.

`exampleData` is a list containing 3 elements:

  - `exampleData$genotypes`: `data.frame` containing the genotypic data
    encoded in allele doses of 100 fictitious individuals with 3333 SNP
    markers. These individuals have 10 chromosomes of length 10^6 bases
    pairs.
  - `exampleData$snpCoord`: `data.frame` containing the coordinates of
    the 3333 individuals’ markers. This data.frame contains 3 columns:
    `chr`, `pos` and `SNPid`.
  - `exampleData$snpEffects`: `numeric` vector containing the “ture”
    effects of the 3333 individuals’ markers about a fictitious
    quantitative trait based on an additive architecture.

#### initialisation

First we must load the package:

``` r
library(breedSimulatR)
#> Loading required package: R6
```

#### Specie specification

Let’s specify the specie:

``` r
# create specie object
specie_statEx <- specie$new(specName = "Statisticae exempli",
                        nChr = 10,
                        lchr = 1e6,
                        ploidy = 2,
                        recombRate = 3/1e6)
#> A new species has emerged: Statisticae exempli !
```

#### SNP specification

We must specify the information about the positions of the genotypic
markers used in the simulation.

Let’s load these information (stored in `exampleData$snpCoord`) and
create the `SNPinfo` object.

``` r
# data preview
head(exampleData$snpCoord)
#>     chr    pos    SNPid
#> 1 Chr10 223757 snp32341
#> 2 Chr01 937638 snp03760
#> 3 Chr05 933373 snp17846
#>  [ reached 'max' / getOption("max.print") -- omitted 3 rows ]
```

``` r
# create SNPinfo object
SNPs <- SNPinfo$new(SNPcoord = exampleData$snpCoord,
                    specie = specie_statEx)
print(SNPs)
#> specie: Statisticae exempli
#> 3333 Markers on 10 chromosomes :
#> Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 
#>   415   247   425   322   381   269   238   355   342   339 
#> SNPcoord:
#>        chr  pos    SNPid
#> 338  Chr01 2068 snp00006
#> 1450 Chr01 2708 snp00009
#> 1391 Chr01 2782 snp00011
#>  [ reached 'max' / getOption("max.print") -- omitted 3330 rows ]
```

#### Population initialisation

We can now generate an initial population from genotypic data.

Let’s load the genotypic information (stored in `exampleData$genotypes`)
and create the `population` object:

``` r
# data preview
exampleData$genotypes[1:3,1:5]
#>          snp00006 snp00009 snp00011 snp00018 snp00026
#> Coll0001        2        2        2        0        2
#> Coll0002        0        2        2        2        0
#>  [ reached 'max' / getOption("max.print") -- omitted 1 rows ]
```

``` r
# create population object
initPop <- createPop(geno = exampleData$genotypes,
                     SNPinfo = SNPs,
                     popName = "Initial population")
```

#### Selection Simulation

In order to perform crossing, we must specify wich individuals must be
mate together. Therefore, we must create functions which generate a
crossing table from our population.

For this example, we will use the function `selectBV`, which returns the
names of the best individuals according to their breeding values.

Then, the function `randomMate` will generate the crossing table.

This function must return a data.frame specifying the crosses.

``` r
exampleData$snpEffects
#>    snp00006    snp00009    snp00011    snp00018    snp00026    snp00031 
#>  0.05110930  0.08099233  0.18353259  0.02819434  0.14413695 -0.04454622 
#>    snp00035    snp00036    snp00049    snp00052 
#> -0.06519311  0.04378479  0.15264373 -0.19802910 
#>  [ reached getOption("max.print") -- omitted 3323 entries ]
(selectedInds <- selectBV(pop = initPop,
                          SNPeffects = exampleData$snpEffects,
                          n = 10))
#>  [1] "Coll0028" "Coll0063" "Coll0034" "Coll0086" "Coll0059" "Coll0080"
#>  [7] "Coll0070" "Coll0021" "Coll0087" "Coll0027"

(crossTable <- randomMate(inds = selectedInds,
                          n = 120,
                          names = "generation_1"))
#>       ind1     ind2 n            names
#> 1 Coll0059 Coll0028 1 generation_1-001
#> 2 Coll0080 Coll0080 1 generation_1-002
#>  [ reached 'max' / getOption("max.print") -- omitted 118 rows ]
```

We can now generate the offsprings:

``` r
newPop <- population$new(name = "1st offsprings",
                         inds = makeCrosses(crosses = crossTable, pop = initPop))
```

``` r
newPop
#> Population: 1st offsprings
#> Species: Statisticae exempli
#> Number of individuals: 120
```

This process can be included in loops in order to simulate several
generations.

## Issues

When encountering a problem with the package or if you have questions,
please report issues on GitHub
[here](https://github.com/ut-biomet/breedSimulatR/issues).

I will do my best to help you as soon as possible.

## Contributing

You can contribute in various ways:

  - report an [issue](https://github.com/ut-biomet/breedSimulatR/issues)
    (online, see the above section)
  - suggest improvements (in the same way as issues)
  - propose a pull request (after creating a new branch)

When editing the content of this package, please run the following
commands before asking a pull request:

``` r
devtools::document()
pkg <- devtools::build()
devtools::check_built(pkg)
```

## Citation

Please cite this package when using it for your projects:

``` r
citation("breedSimulatR")
```

See also `citation()` for citing R itself.

## Acknowledgments

Thanks to [Kosuke Hamazaki](https://github.com/KosukeHamazaki) for his
feedbacks.

## References

`breedSimulatR` is written in **R**:

<ul>

<li>

R Core Team (2019).<em>R: A Language and Environment for Statistical
Computing</em>.R Foundation for Statistical Computing, Vienna,
Austria.<a href="https://www.R-project.org/">https://www.R-project.org/</a>.

</li>

</ul>

`breedSimulatR` package or its devellopment required the following R
packages:

<ul>

<li>

<b> R6 </b>

<ul>

<li>

Chang W (2019).<em>R6: Encapsulated Classes with Reference
Semantics</em>.R package version 2.4.1,
<a href="https://CRAN.R-project.org/package=R6">https://CRAN.R-project.org/package=R6</a>.

</li>

</ul>

</li>

<li>

<b> plotly </b>

<ul>

<li>

Sievert C (2018).<em>plotly for
R</em>.<a href="https://plotly-r.com">https://plotly-r.com</a>.

</li>

</ul>

</li>

<li>

<b> devtools </b>

<ul>

<li>

Wickham H, Hester J, Chang W (2019).<em>devtools: Tools to Make
Developing R Packages Easier</em>.R package version 2.2.1,
<a href="https://CRAN.R-project.org/package=devtools">https://CRAN.R-project.org/package=devtools</a>.

</li>

</ul>

</li>

<li>

<b> roxygen2 </b>

<ul>

<li>

Wickham H, Danenberg P, CsÃ¡rdi G, Eugster M (2019).<em>roxygen2:
In-Line Documentation for R</em>.R package version 7.0.2,
<a href="https://CRAN.R-project.org/package=roxygen2">https://CRAN.R-project.org/package=roxygen2</a>.

</li>

</ul>

</li>

<li>

<b> testthat </b>

<ul>

<li>

Wickham H (2011).“testthat: Get Started with Testing.”<em>The R
Journal</em>, <b>3</b>,
5–10.<a href="https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf">https://journal.r-project.org/archive/2011-1/RJournal\_2011-1\_Wickham.pdf</a>.

</li>

</ul>

</li>

<li>

<b> covr </b>

<ul>

<li>

Hester J (2019).<em>covr: Test Coverage for Packages</em>.R package
version 3.4.0,
<a href="https://CRAN.R-project.org/package=covr">https://CRAN.R-project.org/package=covr</a>.

</li>

</ul>

</li>

<li>

<b> knitr </b>

<ul>

<li>

Xie Y (2020).<em>knitr: A General-Purpose Package for Dynamic Report
Generation in R</em>.R package version 1.27,
<a href="https://yihui.org/knitr/">https://yihui.org/knitr/</a>.

</li>

<li>

Xie Y (2015).<em>Dynamic Documents with R and knitr</em>, 2nd
edition.Chapman and Hall/CRC, Boca Raton, Florida.ISBN 978-1498716963,
<a href="https://yihui.org/knitr/">https://yihui.org/knitr/</a>.

</li>

<li>

Xie Y (2014).“knitr: A Comprehensive Tool for Reproducible Research in
R.”In Stodden V, Leisch F, Peng RD (eds.), <em>Implementing Reproducible
Computational Research</em>.Chapman and Hall/CRC.ISBN 978-1466561595,
<a href="http://www.crcpress.com/product/isbn/9781466561595">http://www.crcpress.com/product/isbn/9781466561595</a>.

</li>

</ul>

</li>

<li>

<b> rmarkdown </b>

<ul>

<li>

Allaire J, Xie Y, McPherson J, Luraschi J, Ushey K, Atkins A, Wickham H,
Cheng J, Chang W, Iannone R (2020).<em>rmarkdown: Dynamic Documents for
R</em>.R package version 2.1,
<a href="https://github.com/rstudio/rmarkdown">https://github.com/rstudio/rmarkdown</a>.

</li>

<li>

Xie Y, Allaire J, Grolemund G (2018).<em>R Markdown: The Definitive
Guide</em>.Chapman and Hall/CRC, Boca Raton, Florida.ISBN 9781138359338,
<a href="https://bookdown.org/yihui/rmarkdown">https://bookdown.org/yihui/rmarkdown</a>.

</li>

</ul>

</li>

</ul>

Example data were generated using the serious game “PlantBreedGame”
available on GitHub <https://github.com/timflutre/PlantBreedGame>

<ul>

<li>

<b>PlantBreedGame</b>

<ul>

<li>

Flutre T, Diot J, and David J (2019). <em>PlantBreedGame</em>: A Serious
Game that Puts Students in the Breeder’s Seat. <em>Crop Science.</em>
[DOI 10.2135/cropsci2019.03.0183le](https://dl.sciencesocieties.org/publications/cs/abstracts/59/4/1374)

</li>

</ul>

</li>

</ul>

## License and copyright :copyright:

The `breedSimulatR` package as a whole is licensed under the MIT. See
the [LICENSE](LICENSE) file for more details.

:copyright: The copyright holder is [The Univerity of
Tokyo](https://www.u-tokyo.ac.jp/en/), Laboratory of Biometry and
Bioinformatics.

<img src="man/figures/UTlogo-sm.jpg" style="display: block; margin: auto;"/>
