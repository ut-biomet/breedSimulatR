# breedSimulatR 0.3.0

## Breaking changes

- For the initialization of the class `specie`, parameters `ploidy` and `mutRate` had been removed.
- For functions `selectBV` and `selectWBV`, parameter `SNPeffect` have been renamed to `QTNeffect`.

## Improved Documentation

- Create an [**example** vignette](https://ut-biomet.github.io/breedSimulatR/articles/example.html) presenting a use case example for the package.

# breedSimulatR 0.2.0

## Improvements

- The initialization of the object `specie` now take as input chromosomes length in centimorgans: `lchrCm`. The argument `recombRate` is now deprecated.
- Crossing over simulations can now take in account the linkage map position of the markers. 

## Bug Fixes

- `specie`'s `getChrLength` method now works fine with input of class "character"


# breedSimulatR 0.1.3

## Improvements

- The method `phenotyper$trial` can now manage the parameter `rep` as a vector for each individual.

# breedSimulatR 0.1.2

## Feature

- Phenotyping features implemented with the objects `trait` and `phenotyper`.

# breedSimulatR 0.1.1

## Bug Fixes

- Better unique name generation for descendants

# breedSimulatR 0.1.0

First version of the breedSimulatR package.
