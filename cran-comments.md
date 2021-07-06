
## Resubmission

This is a resubmission.

Since the last submission (published 2021-02-04), I changed the following:

* The function `wlsPower()` now has an argument `alpha_012` that offers an alternative
way to specifiy the correlation matrix.
* In function `wlsPower()`, the argument `AR` now accepts a vector of up to three values. 
This allows to specifiy autoregressive structures for only a subset of: random cluster intercept, 
random intervention effect and random subject intercept. 
* Closed formulae were added. 
* The method `plot.wlsPower` now produces up to three plots, the hatmatrix, the intervention design and the covariance matrix.
* The vignette was extended.


## Test environments

* local R installation , R 4.1.0
* Mac OS 10.15.7 (on GitHub Actions), 4.1.0
* ubuntu 20.04.2 (on GitHub Actions), 4.1.0
* Windows Server 2019 10.0.17763 (on GitHub Actions), R 4.1.0

## R CMD check results

0 errors | 0 warnings | 0 notes
