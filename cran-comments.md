
## Submission of version 0.4.0

### Changes to previous version:

* `glsPower()` now supports count outcomes via `family="poisson"`
* Added ICC (intracluster correlation) transformation functions: `icc_to_RandEff()`, `RandEff_to_icc()`, `RandEff_to_alpha012()`, and `alpha012_to_RandEff()` for converting between random effects variances and ICC/CAC/IAC parameters
* Added vignette on binomial and count outcomes, with pre-calculated contour plots
* In `glsPower()`, the argument `N` now overrides the `N` stored in a supplied
`DesMat` object
* `glsPower()` now fails gracefully (with a warning) if the information content
cannot be calculated
* Character input arguments (e.g. `dsntype`, `family`) now throw an error if no
known option is sufficiently similar, instead of silently choosing the closest match
* Diagnostic output now uses `message()` instead of `print()`
* Vignette plots now use the plotly partial bundle to reduce package size
* Changed covariance matrix construction to use `fbdiag` (fast block diagonal matrix)
* `plot_CellWeights()` now treats `NA` entries in `incompMat` as unobserved
cluster periods
* Replaced `\()` with `function()` for backward compatibility with older R versions
* Fixed typo in `RandEff_to_alpha`
* Added tests for `construct_DesMat()`, `construct_CovMat()`, and `glsPower()`
* Updated vignettes and improved documentation with additional links in help files
* Roxygen documentation now uses markdown format; re-roxygenised all documentation

### Test environments

* local R installation (Linux Mint 22.3), R 4.6.1
* Mac OS (on GitHub Actions), R release
* Windows Server (on GitHub Actions), R release
* ubuntu 22.04 (on GitHub Actions), R devel, release and oldrel-1

### R CMD check results

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
```


## Submission of version 0.3.5

### Changes to previous version:

* internal restructuring and error fixes
* `plot.DesMat()` now has a new option to plot the treatment allocation for individuals 
(instead of entire clusters)
 
 
### Test environments 

* local R installation (Win 10), R 4.4.0
* Mac OS 14.4.1 (on GitHub Actions), R 4.4.0

### R CMD check results

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
```


## Submission of version 0.3.4 

* Fixed the roxygen2 bug as explained in https://github.com/r-lib/roxygen2/issues/1491 


## Submission of version 0.3.2

### Changes to version 0.3.1 

* The most noticeable change in this version is that the abbrevation `wls` 
(weighted least squares) in function names is now replaced with `gls`
(generalised least squares) to more properly reflect the scope of the functionality.
For example, the function `wlsPower()` is now called `glsPower()` - although the
former version still works and throws a warning. 
* The closed formula for the computation of information content is now a dedicated formula, 
called `compute_InfoContent()` 
* In `plot.glsPower()` there now is an option to manually set the font size of the
annotation in the influence plots


### Test environments

* local R installation (Windows 10) , R-devel (4.3.0)
* Mac OS 11.6.6 (on GitHub Actions), R 4.2.0
* ubuntu 20.04.4 (on GitHub Actions), R 4.2.0
* Windows Server 2022 10.0.20348 (on GitHub Actions), R 4.2.0

### R CMD check results

```
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
```



## Resubmission of version 0.3.1

You saw:
```
 Found the following (possibly) invalid URLs:
    URL: https://www.unimedizin-mainz.de/imbei/imbei/welcome-page (moved to https://www.unimedizin-mainz.de/imbei/imbei/welcome-page/)
      From: inst/doc/Getting_Started.html 
```

I fixed the trailing slashes in the vignette. 

### Changes in 0.3.1 to previous version

Since submission of version 0.2.0 (published 2021-07-07), I changed the following:

* The function `wlsPower()` now also computes the information content of 
cluster-period cells. Computation is currently done twice, once with a general formula
and once explicitly. Information content of whole periods or clusters is also computed.
* The method `plot.wlsPower()` recieved multiple updates:
  * It now produces up to four plots: the projection matrix, 
  the information content, the intervention design and the covariance matrix.
  * Incomplete designs (SWD where some cluster-period cells are omitted) are now visualised
  * Plots of projection matrix and information content can now be annotated with particular values in each cell;
  This is the default for smaller designs and can be turned on/off via `annotations = <TRUE/FALSE>`
  * An option `show_colorbar` to hide colour bars was added
  * An option `marginal_plots` to hide marginal plots on whole periods or clusters was added.
  * Various aesthetic improvements, e.g.: Improved hover information, dynamic gap size between cells.
* Vignette was extended

### Test environments

* local R installation (Windows 10) , R 4.1.2
* Mac OS 10.15.7 (on GitHub Actions), R 4.1.2
* ubuntu 20.04.2 (on GitHub Actions), R 4.1.2
* Windows Server 2019 10.0.17763 (on GitHub Actions), R 4.1.2

### R CMD check results

0 errors | 0 warnings | 1 note

```
checking installed package size ... NOTE
  installed size is  5.1Mb
  sub-directories of 1Mb or more:
    doc   5.0Mb
```
The main culprit are the interactive plots, produced with `plotly`. If necessary, I
could remove those from the vignette, but I'd be thrilled if the vignette could keep those plots. 
