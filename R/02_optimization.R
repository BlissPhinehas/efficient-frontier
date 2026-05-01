# 02_optimization.R
# Purpose: Compute expected returns and covariance matrix,
#          run quadratic optimization to trace the efficient frontier,
#          and compute the Capital Market Line (CML).

library(quadprog)
library(zoo)

# ── 1. Load data ──────────────────────────────────────────────────────────────

log_returns <- readRDS("outputs/data/returns_full.rds")

# ── 2. Annualized expected returns and covariance matrix ──────────────────────

mu    <- colMeans(log_returns) * 252          # annualized expected returns
Sigma <- cov(log_returns) * 252               # annualized covariance matrix
n     <- length(mu)                           # number of assets

cat("Annualized Expected Returns:\n")
print(round(mu, 4))
cat("\nCovariance Matrix (first 4x4):\n")
print(round(Sigma[1:4, 1:4], 6))

# ── 3. Quadratic optimization — trace the efficient frontier ──────────────────
# We solve the standard Markowitz QP:
#   minimize    (1/2) w' Sigma w
#   subject to  mu' w  = target_return
#               sum(w) = 1
#               w >= 0   (long-only constraint)
#
# quadprog::solve.QP expects:
#   Dmat = covariance matrix (what we're minimizing)
#   dvec = zero vector
#   Amat = constraint matrix (transposed)
#   bvec = constraint RHS
#   meq  = number of equality constraints (first meq rows of Amat)

optimize_portfolio <- function(target_return, mu, Sigma, long_only = TRUE) {
  n    <- length(mu)
  Dmat <- 2 * Sigma                          # quadprog expects 2*Sigma
  dvec <- rep(0, n)

  # Equality constraints: return target + weights sum to 1
  A_eq <- rbind(mu, rep(1, n))              # 2 x n
  b_eq <- c(target_return, 1)

  if (long_only) {
    # Inequality constraints: w_i >= 0
    A_ineq <- diag(n)                        # n x n identity
    Amat   <- t(rbind(A_eq, A_ineq))        # quadprog wants transposed
    bvec   <- c(b_eq, rep(0, n))
    meq    <- 2
  } else {
    Amat <- t(A_eq)
    bvec <- b_eq
    meq  <- 2
  }

  result <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = meq),
    error = function(e) NULL
  )

  if (is.null(result)) return(NULL)

  weights <- result$solution
  port_return <- sum(weights * mu)
  port_vol    <- sqrt(t(weights) %*% Sigma %*% weights)

  list(
    weights = weights,
    return  = port_return,
    vol     = as.numeric(port_vol)
  )
}

# ── 4. Sweep target returns across feasible range ─────────────────────────────

target_returns <- seq(min(mu) * 0.8, max(mu) * 1.0, length.out = 200)

frontier <- data.frame(vol = numeric(), ret = numeric())

for (tr in target_returns) {
  port <- optimize_portfolio(tr, mu, Sigma)
  if (!is.null(port)) {
    frontier <- rbind(frontier, data.frame(vol = port$vol, ret = port$ret))
  }
}

# Keep only the upper half (efficient frontier, not the dominated lower half)
min_vol_idx    <- which.min(frontier$vol)
efficient_frontier <- frontier[min_vol_idx:nrow(frontier), ]

cat("\nFrontier points computed:", nrow(frontier), "\n")
cat("Efficient frontier points:", nrow(efficient_frontier), "\n")

# ── 5. Minimum variance portfolio ─────────────────────────────────────────────

mvp <- frontier[min_vol_idx, ]
cat("\nMinimum Variance Portfolio:\n")
cat("  Volatility:", round(mvp$vol, 4), "\n")
cat("  Return:    ", round(mvp$ret, 4), "\n")

# ── 6. Capital Market Line ─────────────────────────────────────────────────────
# Risk-free rate: approximate 10-year average Fed Funds rate for 2000-2023

rf <- 0.02  # 2% annualized risk-free rate

# Tangency portfolio = max Sharpe ratio portfolio
sharpe_ratios <- (frontier$ret - rf) / frontier$vol
tangency_idx  <- which.max(sharpe_ratios)
tangency      <- frontier[tangency_idx, ]

cat("\nTangency Portfolio (Max Sharpe):\n")
cat("  Volatility:", round(tangency$vol, 4), "\n")
cat("  Return:    ", round(tangency$ret, 4), "\n")
cat("  Sharpe:    ", round((tangency$ret - rf) / tangency$vol, 4), "\n")

# CML: line from (0, rf) through tangency portfolio
cml_vols <- seq(0, max(frontier$vol) * 1.3, length.out = 100)
cml_rets <- rf + ((tangency$ret - rf) / tangency$vol) * cml_vols
cml       <- data.frame(vol = cml_vols, ret = cml_rets)

# ── 7. Individual asset risk/return ───────────────────────────────────────────

asset_points <- data.frame(
  ticker = colnames(log_returns),
  vol    = sqrt(diag(Sigma)),
  ret    = mu
)

cat("\nIndividual Asset Risk/Return:\n")
print(round(asset_points[, c("vol", "ret")], 4))

# ── 8. Save results ───────────────────────────────────────────────────────────

saveRDS(frontier,           "outputs/data/frontier_full.rds")
saveRDS(efficient_frontier, "outputs/data/efficient_frontier.rds")
saveRDS(cml,                "outputs/data/cml.rds")
saveRDS(tangency,           "outputs/data/tangency.rds")
saveRDS(asset_points,       "outputs/data/asset_points.rds")
saveRDS(list(mu = mu, Sigma = Sigma, rf = rf), "outputs/data/params.rds")

cat("\nOptimization complete. Results saved to outputs/data/\n")