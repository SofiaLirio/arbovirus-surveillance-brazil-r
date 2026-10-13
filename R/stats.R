# Statistical routines written to reproduce the Python libraries the project was first built with:
#   mk_original_test()  pymannkendall.original_test (+ sens_slope)
#   spearman_test()     scipy.stats.spearmanr (two-sided, t approximation)
#   moran_global()      esda.Moran (row-standardised weights, 999 permutations, folded p_sim)
#   moran_local()       esda.Moran_Local (conditional permutation, folded p_sim, quadrants 1-4)
# Analytic quantities are identical to the Python ones; permutation p-values use R's RNG, so they
# agree with the Python values only up to Monte-Carlo error.

# ---- Mann-Kendall ------------------------------------------------------------------------------
mk_original_test <- function(x_old, alpha = 0.05) {
  x_old <- as.numeric(x_old)
  x <- x_old[!is.na(x_old)]
  n <- length(x)
  s <- 0
  if (n > 1) for (k in seq_len(n - 1)) s <- s + sum(x[(k + 1):n] > x[k]) - sum(x[(k + 1):n] < x[k])
  ux <- unique(x)
  if (length(ux) == n) {
    var_s <- (n * (n - 1) * (2 * n + 5)) / 18
  } else {
    tp <- as.numeric(table(match(x, ux)))
    var_s <- (n * (n - 1) * (2 * n + 5) - sum(tp * (tp - 1) * (2 * tp + 5))) / 18
  }
  tau <- s / (0.5 * n * (n - 1))
  z <- if (s > 0) (s - 1) / sqrt(var_s) else if (s == 0) 0 else (s + 1) / sqrt(var_s)
  p <- 2 * (1 - stats::pnorm(abs(z)))
  h <- abs(z) > stats::qnorm(1 - alpha / 2)
  trend <- if (z < 0 && h) "decreasing" else if (z > 0 && h) "increasing" else "no trend"
  ss <- sens_slope(x_old)
  list(trend = trend, h = h, p = p, z = z, Tau = tau, s = s, var_s = var_s,
       slope = ss$slope, intercept = ss$intercept)
}

# Theil-Sen slope over index steps (NaN kept in place) and Conover intercept anchored at index 0
sens_slope <- function(x) {
  x <- as.numeric(x); n <- length(x)
  d <- numeric(0)
  if (n > 1) for (i in seq_len(n - 1)) {
    j <- (i + 1):n
    d <- c(d, (x[j] - x[i]) / (j - i))
  }
  slope <- if (all(is.na(d))) NA_real_ else stats::median(d, na.rm = TRUE)
  idx <- (seq_len(n) - 1)[!is.na(x)]
  intercept <- stats::median(x, na.rm = TRUE) - stats::median(idx) * slope
  list(slope = slope, intercept = intercept)
}

# ---- Spearman ----------------------------------------------------------------------------------
spearman_test <- function(x, y) {
  n <- length(x)
  rs <- suppressWarnings(stats::cor(rank(x), rank(y)))
  if (is.na(rs)) return(list(rho = NA_real_, p = NA_real_))
  dof <- n - 2
  t <- rs * sqrt(pmax(dof / ((rs + 1) * (1 - rs)), 0))
  p <- 2 * stats::pt(abs(t), dof, lower.tail = FALSE)
  list(rho = rs, p = p)
}

# ---- Natural breaks ----------------------------------------------------------------------------
# Jenks upper class bounds (data values): the optimum mapclassify.NaturalBreaks / FisherJenks return.
natural_breaks <- function(vals, k) {
  k <- min(k, length(unique(vals)))
  b <- classInt::classIntervals(vals, n = k, style = "jenks")$brks
  list(bins = b[-1], k = length(b) - 1, min = b[1])
}
# np.digitize(vals, bins, right=True): number of bins strictly below each value
digitize_right <- function(vals, bins) vapply(vals, function(v) sum(bins < v), integer(1))

# ---- Spatial weights ---------------------------------------------------------------------------
# Row-standardised KNN weights as a list: nb[[i]] neighbour indices, wt[[i]] weights.
knn_listw <- function(coords, k = 5) {
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k = k, longlat = FALSE))
  list(nb = lapply(nb, as.integer), wt = lapply(nb, function(v) rep(1 / length(v), length(v))), n = length(nb),
       spdep = spdep::nb2listw(nb, style = "W"))
}
spatial_lag <- function(w, z) vapply(seq_len(w$n), function(i) sum(w$wt[[i]] * z[w$nb[[i]]]), numeric(1))

folded_p <- function(sim, obs, perms) {
  larger <- sum(sim >= obs)
  if ((perms - larger) < larger) larger <- perms - larger
  (larger + 1) / (perms + 1)
}

# ---- Global Moran's I --------------------------------------------------------------------------
moran_global <- function(y, w, permutations = 999, seed = NULL) {
  y <- as.numeric(y); n <- length(y)
  z <- y - mean(y)
  z2ss <- sum(z * z)
  s0 <- sum(unlist(w$wt))
  calc <- function(zz) n / s0 * sum(zz * spatial_lag(w, zz)) / z2ss
  I <- calc(z)
  out <- list(I = I, EI = -1 / (n - 1))
  if (permutations) {
    if (!is.null(seed)) set.seed(seed)
    sim <- vapply(seq_len(permutations), function(i) calc(sample(z)), numeric(1))
    sd0 <- sqrt(mean((sim - mean(sim))^2))
    out$p_sim <- folded_p(sim, I, permutations)
    out$EI_sim <- mean(sim)
    out$z_sim <- (I - mean(sim)) / sd0
    out$sim <- sim
  }
  out
}

# ---- Local Moran's I (LISA) --------------------------------------------------------------------
moran_local <- function(y, w, permutations = 999, seed = NULL) {
  y <- as.numeric(y); n <- length(y)
  z <- (y - mean(y)) / sqrt(mean((y - mean(y))^2))
  den <- sum(z * z)
  zl <- spatial_lag(w, z)
  Is <- (n - 1) * z * zl / den
  zp <- z > 0; lp <- zl > 0
  q <- ifelse(zp & lp, 1L, ifelse(!zp & lp, 2L, ifelse(!zp & !lp, 3L, 4L)))
  p_sim <- rep(NA_real_, n)
  if (permutations) {
    if (!is.null(seed)) set.seed(seed)
    scaling <- (n - 1) / den
    for (i in seq_len(n)) {
      others <- z[-i]; k <- length(w$nb[[i]]); wi <- w$wt[[i]]
      draws <- matrix(replicate(permutations, sample(others, k)), nrow = k)
      sim <- z[i] * colSums(draws * wi) * scaling
      p_sim[i] <- folded_p(sim, Is[i], permutations)
    }
  }
  list(Is = Is, q = q, p_sim = p_sim, z = z)
}
