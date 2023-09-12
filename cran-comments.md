## Submission of version 0.3.4 

* Fixed the roxygen2 bug as explained in https://github.com/r-lib/roxygen2/issues/1491 


## Submission of version 0.3.2

### Changes to version 0.3.1 

* The most noticeable change in this version is that the abbrevation `wls` 
(weighted least squares) in function names is now replaced with `gls`
(generalised least squars) to more properly reflect the scope of the functionality.
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
