#' @title Power distribution over cluster-to-sequence allocations
#'
#' @description
#' Wrapper around \code{\link{glsPower}} for designs in which the clusters
#' have *known but unequal* sizes. Because the randomization assigns clusters
#' to sequences, the resulting power is a random variable over the set of
#' possible allocations.
#'
#' `N_alloc_Power` enumerates (or samples) all distinct ways of assigning the
#' supplied clusters to the sequences defined by `Cl`, computes the power for
#' each allocation via \code{\link{glsPower}}, and summarises the resulting
#' distribution (mean, min, max, quartiles, ...).
#'
#' Within a sequence, the order of clusters does not affect power (cluster
#' contributions to the information add up commutatively), so an *allocation*
#' is a partition of the cluster set into the sequences, not a full
#' permutation. Allocations that differ only by swapping equally-sized
#' clusters therefore share the same power; `N_alloc_Power` de-duplicates
#' them and weights accordingly.
#'
#' @param N_pool numeric vector of length `sum(Cl)`, number of
#' individuals per cluster, i.e. the pool of per-cluster `N` values *before*
#' allocation to sequences. The order of this vector defines the cluster
#' identities (cluster 1, 2, ...); the function averages over all ways of
#' assigning these clusters to sequences.
#' @param Cl integer (vector), number of clusters per sequence group,
#' as in \code{\link{glsPower}}. Must satisfy `sum(Cl) == length(N_pool)`.
#' @param ... further arguments passed to \code{\link{glsPower}}, e.g. `mu0`,
#' `mu1`, `sigma`, `tau`, `family`, `alpha_0_1_2`, `eta`, `rho`, `gamma`,
#' `psi`, `AR`, `sig.level`, `trtDelay`, `incomplete`, `timeAdjust`,
#' `dsntype`, `timepoints`, `INDIV_LVL`. Arguments `N` and `power` are
#' managed by `N_alloc_Power` and must NOT be supplied via `...`.
#' @param n_MC integer, number of Monte Carlo allocations to draw uniformly at
#' random. If `NULL` (the default), all distinct allocations are enumerated
#' exhaustively, subject to `max_alloc`. Use this when the number of distinct
#' allocations is too large to enumerate.
#' @param max_alloc integer, upper bound on the number of *unique*
#' allocations evaluated in exhaustive mode. If the number of unique
#' allocations exceeds this, the function stops with an error suggesting
#' to set `n_MC`. The underlying enumeration of index-partitions is
#' capped internally at `1e6`.
#' @param parallel logical, should the per-allocation \code{\link{glsPower}}
#' calls be evaluated in parallel via \code{future.apply::future_mapply}?
#' If `TRUE`, package `future.apply` is required and the user must select
#' a parallel backend beforehand, e.g.
#' `future::plan(future::multisession)`. If `FALSE` (the default), a
#' plain `mapply` is used.
#' @param keep character, one of `"summary"` (the default) or `"all"`. The
#' former returns only the summary statistics, the latter additionally the
#' per-allocation powers and the corresponding cluster-size assignments.
#' @param verbose logical, print progress messages.
#'
#' @details
#' Let \eqn{N = (N_1, \dots, N_K)} the vector of cluster sizes and
#' \eqn{\theta:= \mu_1-\mu_0} the treatment effect under investigation.
#' For each allocation of the clusters to the sequences, the power of the
#' Wald test of \eqn{\theta} is computed by \code{\link{glsPower}}.
#' The function returns summary statistics of the resulting power
#' distribution over all (equally weighted) allocations.
#'
#' @return
#' An object of class `N_alloc_Power` (a list). The `summary` element is a
#' named numeric vector with `n_allocations`, `n_unique`, `mean`, `sd`,
#' `min`, `Q1`, `median`, `Q3`, `max` and `IQR`. The `power` element gives the
#' (unique) power values and `weights` their multiplicities in exhaustive
#' mode (all `1` in Monte Carlo mode, i.e. equal weights); these two together
#' define the power distribution and can be plotted with
#' \code{\link[=plot.N_alloc_Power]{plot}}. The `best` and `worst` elements
#' give the cluster-size assignment (N-vector, ordered by sequence) achieving
#' the maximum and minimum power. If `keep="all"`, `N_allocations`
#' additionally contains the unique cluster-size assignments as rows of a
#' matrix.
#'
#' @export
#'
#' @examples
#' ## 8 clusters of unequal size, 4 sequences with 2 clusters each.
#' ## Exhaustive enumeration (2520 index-partitions).
#' N_alloc_Power(N_pool = c(20, 30, 25, 15, 40, 10, 35, 22),
#'               Cl = c(2, 2, 2, 2),
#'               mu0 = 0.2, mu1 = 0.35, tau = 0.5, family = "binomial")
#'
#' ## Gaussian outcome, vary cluster sizes
#' N_alloc_Power(N_pool = c(12, 8, 10, 9, 14),
#'               Cl = rep(1, 5),
#'               mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33)
#'
#' ## Larger design: use Monte Carlo when enumeration is infeasible.
#' \dontrun{
#' set.seed(123)
#' future::plan(future::multisession)
#' N_alloc_Power(N_pool = round(runif(16, 10, 40)),
#'               Cl = rep(2, 8),
#'               mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33,
#'               n_MC = 5000, parallel = TRUE)
#' }
#'
N_alloc_Power <- function(N_pool, Cl, ...,
                        n_MC = NULL, max_alloc = 50000,
                        parallel = FALSE,
                        keep = "summary",
                        verbose = TRUE) {

  keep   <- choose_character_Input(c("summary", "all"), keep)
  Cl     <- as.integer(Cl)
  N_pool <- as.numeric(N_pool)
  sumCl  <- length(N_pool)
  lenCl  <- length(Cl)

  if (sum(Cl) != sumCl)
    stop("sum(Cl) (=", sum(Cl), ") must equal length(N_pool) (=", sumCl, ").")
  if (any(Cl <= 0))     stop("All entries of Cl must be positive.")
  if (any(N_pool <= 0)) stop("All entries of N_pool must be positive.")

  ## Check additional args passed to glsPower
  dots <- list(...)
  if ("N"     %in% names(dots)) stop("'N' is derived from 'N_pool' and managed by N_alloc_Power().")
  if ("power" %in% names(dots)) stop("'power' is a random variable in N_alloc_Power().")

  #### Build the set of (N-vector, weight) to evaluate. ####
  if (is.null(n_MC)) {
    ## Exhaustive enumeration of all index-partitions. De-duplicate the
    ## resulting N-vectors (equal cluster sizes) and weight by multiplicity.
    multinom_coeff <- round( exp(lgamma(sumCl + 1) - sum(lgamma(Cl + 1))) ) 
    ## raw enumeration cap: the label matrix and N_all must fit into memory
    if (multinom_coeff > 1e6)
      stop(multinom_coeff, " index-partitions to enumerate. Use n_MC for Monte Carlo sampling instead.")
    if (verbose)
      message("Enumerating allocations (~", round(multinom_coeff), " index-partitions)...")

    ## All distinct assignments of clusters to sequences: the permutations
    ## of the multiset rep(1:lenCl, Cl), each appearing exactly once.
    ## N ordered by sequence, sizes sorted within (for the deduplication).
    labs    <- RcppAlgos::permuteGeneral(lenCl, sumCl, freqs = Cl)
    N_all   <- lapply(seq_len(nrow(labs)), function(r) N_pool[order(labs[r, ], N_pool)])
    keys    <- vapply(N_all, paste, collapse = ",", character(1))
    dup     <- collapse::fduplicated(keys)
    weights <- as.integer(collapse::qtab(keys)[keys[!dup]])
    N_list  <- N_all[!dup]
    ## max_alloc guards the expensive part: the per-allocation glsPower calls
    if (length(N_list) > max_alloc)
      stop(length(N_list), "unique allocations to evaluate, which exceeds 'max_alloc' (=", max_alloc, ").\n",
           "Set n_MC to a positive integer to use Monte Carlo sampling ",
           "instead, or increase max_alloc (at your own risk).")
    mode    <- "exhaustive"
    if (verbose)  message("Found ", length(N_list), " unique allocation(s).")
  } else {
    ## Monte Carlo: draw n_MC allocations uniformly at random, store them in N_list
    n_MC <- as.integer(n_MC)
    if (!isTRUE(n_MC > 0L)) stop("n_MC must be a positive integer.")
    N_list  <- replicate(n_MC, N_pool[sample(1:sumCl)], simplify = FALSE)
    weights <- rep(1, n_MC)
    mode    <- "monte_carlo"
    if (verbose) message("Sampling ", n_MC, " Monte Carlo allocation(s)...")
  }

  ## Get power for each N-vector in N_list via glsPower. `MoreArgs` holds the constant arguments.
  if (verbose) message("Computing power for ", length(N_list), " allocation(s)...")

  MoreArgs <- c(list(Cl = Cl, verbose = 0), dots)
  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("Package 'future.apply' is required for parallel=TRUE. ",
           "Install it and set a backend, e.g. future::plan(future::multisession).")
    powers <- future.apply::future_mapply(glsPower, N = N_list,
                                          MoreArgs = MoreArgs)
  } else {
    powers <- mapply(glsPower, N = N_list, MoreArgs = MoreArgs)
  }

  #### Summary statistics -----####
  qs <- collapse::fquantile(powers, w = weights, names = FALSE)
  summ <- c(n_allocations = sum(weights),
            n_unique      = if (mode == "exhaustive") length(N_list) else NA_integer_,
            mean = collapse::fmean(powers, w = weights),
            sd   = collapse::fsd(  powers, w = weights),
            min  = qs[1], Q1  = qs[2], median = qs[3],
            Q3   = qs[4], max = qs[5], IQR    = qs[4] - qs[2])

  out <- list(summary = summ, power = powers, weights = weights,
              N_allocations = if (keep == "all")
                matrix(unlist(N_list), nrow = length(N_list), byrow = TRUE) else NULL,
              best  = N_list[[which.max(powers)]],
              worst = N_list[[which.min(powers)]],
              mode  = mode, n_MC = n_MC,
              Cl    = Cl  , N_pool = N_pool)
  class(out) <- "N_alloc_Power"
  out
}


#' @title Print an object of class `N_alloc_Power`
#'
#' @param x object of class N_alloc_Power
#' @param ... Arguments to be passed to methods
#'
#' @method print N_alloc_Power
#'
#' @return Invisibly returns `x`.
#'
#' @export
#'
print.N_alloc_Power <- function(x, ...) {
  cat("<N_alloc_Power> power distribution over cluster allocations\n")
  cat("  mode            :", x$mode, "\n")
  if (x$mode == "monte_carlo")
    cat("  Monte Carlo draws:", x$n_MC, "\n")
  cat("  Cl              :", paste(x$Cl, collapse = ", "), "\n")
  cat("  N_pool          :", paste(x$N_pool, collapse = ", "), "\n")
  cat("  allocations     :", x$summary[["n_allocations"]], "\n")
  if (!is.na(x$summary[["n_unique"]]))
    cat("  unique allocs   :", x$summary[["n_unique"]], "\n")
  cat("\n  Power summary:\n")
  print(round(x$summary[c("mean","sd","min","Q1","median","Q3","max","IQR")], 4))
  cat("\n  best  (max power) N-vector:", paste(x$best,  collapse = ", "), "\n")
  cat("  worst (min power) N-vector:", paste(x$worst, collapse = ", "), "\n")
  invisible(x)
}


#' @title plot.N_alloc_Power
#'
#' @inheritParams print.N_alloc_Power
#' @param x An object of class `N_alloc_Power`
#' @param bins integer, number of histogram bins, spanning the observed
#' power range plus a small padding.
#' @param show_mean logical, should a dashed vertical line at the mean power
#' be added?
#' @param show_IQR logical, should a shaded band over the interquartile range
#' (Q1 to Q3) be added behind the bars?
#' @param fill character, colour of the histogram bars. Defaults to
#' "steelblue", matching the package's colour scheme.
#'
#' @method plot N_alloc_Power
#'
#' @return a plotly html widget, displaying the power distribution as a
#' histogram with reference lines for the mean and the interquartile range.
#'
#' @export
#'
#' @examples
#' ap <- N_alloc_Power(N_pool = c(12, 8, 10, 9, 14),
#'                     Cl = rep(1, 5),
#'                     mu0 = 0, mu1 = 1, sigma = 2, tau = 0.33)
#' plot(ap)
#'
plot.N_alloc_Power <- function(x, bins = 30, show_mean = TRUE,
                               show_IQR = TRUE, fill = "steelblue", ...){

  powers  <- x$power
  weights <- x$weights

  s      <- x$summary
  rng    <- collapse::frange(powers)
  pad    <- 0.02 * diff(rng)    ## a little horizontal padding so the extreme bars are not clipped
  xrng   <- c(rng[1] - pad, rng[2] + pad)

  shapes <- list()
  ## IQR band, drawn below the bars
  if (isTRUE(show_IQR)) {
    shapes <- c(shapes, list(list(
      type = "rect", xref = "x", yref = "paper",
      x0 = s[["Q1"]], x1 = s[["Q3"]], y0 = 0, y1 = 1,
      fillcolor = "lightgray", line = list(width = 0),
      opacity = 0.45, layer = "below")))
  }
  ## mean, dashed vertical line spanning the plot
  if (isTRUE(show_mean)) {
    shapes <- c(shapes, list(list(
      type = "line", xref = "x", yref = "paper",
      x0 = s[["mean"]], x1 = s[["mean"]], y0 = 0, y1 = 1,
      line = list(color = "firebrick", width = 2, dash = "dash"))))
  }

  annot <- list()
  if (isTRUE(show_mean)) {
    annot <- c(annot, list(list(
      x = s[["mean"]], y = 1, yref = "paper", yanchor = "top",
      text = sprintf("mean = %.3f", s[["mean"]]),
      showarrow = FALSE, font = list(color = "firebrick", size = 12))))
  }

  subtitle <- paste(x$mode, "|", "min", round(s["min"],3) , 
                    "median", round(s["median"],3),  "max", round(s["max"],3))  

  out <- plot_ly(x = ~rep(powers, weights), type = "histogram",
                 xbins = list(start = xrng[1], end = xrng[2],
                             size = diff(xrng) / bins),
                 marker = list(color = fill,
                               line = list(color = "white", width = 0.4))) %>%
    layout(title = list(text = paste0("Power distribution over cluster allocations<br>",
                                      "<sup>", subtitle, "</sup>")),
           xaxis = list(title = "Power", range = xrng),
           yaxis = list(title = "Number of allocations"),
           shapes = shapes,
           annotations = annot,
           bargap = 0.02)
  out
}
