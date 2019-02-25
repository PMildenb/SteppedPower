
# `SteppedPower` - sample size calculation in mixed model settings with focus on stepped wedge designs

`SteppedPower` offers options for power calculation as well as sample size calculation (by the use of uniroot() ) for
longitudinal mixed model settings. All calculations are oracle estimates i.e. assume random effect variances to be known (or guessed) in advance.
A stepped wedge design is thought as a mixed model:
$$y_{ij}= \mu + \alpha_i + X (\theta_{ij} + c_j) + b_j + e_{ij}$$
where  
* $y_{ij}$ is the mean response in cluster $j$ at time $i$
* $\mu$ is a grand total mean
* $\alpha_i$ is the time trend at time $i$
* $X$ is a design matrix (see options below)
The workhorse function of this package is `wlsMixedPower()`.

## Installation
not on CRAN yet, so installation currently via  
`devtools::install_github("PMildenb/SteppedPower")`  
`library(SteppedPower)`

## Design Matrix Options
Currently there are the following, to be specified by the argument `design=`   
* [SWD] stepped wedge designs
* [parallel] parallel designs
* [parallel_baseline] parallel designs with baseline period(s)


## Outcome Options
Currently implemented are
* normal 
* binomial
outcomes


