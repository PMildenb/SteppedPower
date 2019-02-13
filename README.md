
# `SteppedPower` - sample size calculation in mixed model settings with focus on stepped wedge designs

`SteppedPower` offers options for power calculation as well as sample size calculation (by the use of uniroot() ) for
longitudinal mixed model settings. All calculations are oracle estimates i.e. assume random effect variances to be known (or guessed) in advance.
  
The workhorse function of this package is `wlsMixedPower()`.

## Installation
not on CRAN yet, so installation currently via  
`devtools::install_github("PMildenb/SteppedPower")`  
`library(SteppedPower)`

## Design Options
Currently there are the following, to be specified by the argument `design=`  
* [SWD] stepped wedge designs
* [parallel] parallel designs
* [parallel_baseline] parallel designs with baseline period(s)


## Outcome Options
Currently implemented are
* normal 
* binomial
outcomes
