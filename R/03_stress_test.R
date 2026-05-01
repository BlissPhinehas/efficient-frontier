# 03_stress_test.R
# Purpose: Re-run frontier optimization on 2008 and 2020 crisis windows
#          and produce a side-by-side comparison plot showing frontier shift.

library(quadprog)
library(ggplot2)
library(zoo)

# ── 1. Load data ──────────────────────────────────────────────────────────────

returns_full <- readRDS("outputs/data/returns_full.rds")
returns_2008 <- readRDS("outputs/data/returns_2008.rds")
returns_2020 <- readRDS("outputs/data/returns_2020.rds")
params       <- readRDS("outputs/data/params.rds")
rf           <- params$rf

# ── 2. Helper: compute frontier from a returns matrix ─────────────────────────

compute_frontier <- function(returns, label, n_points = 150) {

  mu    <- colMeans(returns) * 252
  Sigma <- cov(returns) * 252
  n     <- length(mu)

  optimize_portfolio <- function(target_return) {
    Dmat   <- 2 * Sigma
    dvec   <- rep(0, n)
    A_eq   <- rbind(mu, rep(1, n))
    b_eq   <- c(target_return, 1)
    A_ineq <- diag(n)
    Amat   <- t(rbind(A_eq, A_ineq))
    bvec   <- c(b_eq, rep(0, n))

    result <- tryCatch(
      solve.QP(Dmat, dvec, Amat, bvec, meq = 2),
      error = function(e) NULL
    )
    if (is.null(result)) return(NULL)
    w      <- result$solution
    list(vol = as.numeric(sqrt(t(w) %*% Sigma %*% w)),
         ret = sum(w * mu))
  }

  targets  <- seq(min(mu) * 0.8, max(mu) * 1.0, length.out = n_points)
  frontier <- do.call(rbind, lapply(targets, function(tr) {
    p <- optimize_portfolio(tr)
    if (!is.null(p)) data.frame(vol = p$vol, ret = p$ret, period = label)
    else NULL
  }))

  # Keep only efficient (upper) half
  min_idx  <- which.min(frontier$vol)
  frontier <- frontier[min_idx:nrow(frontier), ]

  # Compute tangency
  sharpes     <- (frontier$ret - rf) / frontier$vol
  tang_idx    <- which.max(sharpes)
  tangency    <- frontier[tang_idx, ]

  list(frontier = frontier, tangency = tangency,
       mu = mu, Sigma = Sigma, label = label)
}

# ── 3. Compute frontiers for all three periods ────────────────────────────────

cat("Computing full-period frontier...\n")
result_full <- compute_frontier(returns_full, "Full Period (2004–2023)")

cat("Computing 2008 crisis frontier...\n")
result_2008 <- compute_frontier(returns_2008, "2008 Financial Crisis")

cat("Computing 2020 crisis frontier...\n")
result_2020 <- compute_frontier(returns_2020, "2020 COVID Crash")

# ── 4. Print summary stats ────────────────────────────────────────────────────

summarize_result <- function(res) {
  tang <- res$tangency
  mvp  <- res$frontier[which.min(res$frontier$vol), ]
  cat("  Period:          ", res$label, "\n")
  cat("  MVP  vol/ret:    ", round(mvp$vol, 4), "/", round(mvp$ret, 4), "\n")
  cat("  Tang vol/ret:    ", round(tang$vol, 4), "/", round(tang$ret, 4), "\n")
  cat("  Tang Sharpe:     ", round((tang$ret - rf) / tang$vol, 4), "\n\n")
}

cat("\n── Frontier Summary ──────────────────────────────────────────\n\n")
summarize_result(result_full)
summarize_result(result_2008)
summarize_result(result_2020)

# ── 5. Combined plot ──────────────────────────────────────────────────────────

all_frontiers <- rbind(
  result_full$frontier,
  result_2008$frontier,
  result_2020$frontier
)

tangency_points <- rbind(
  data.frame(result_full$tangency, period = "Full Period (2004–2023)"),
  data.frame(result_2008$tangency, period = "2008 Financial Crisis"),
  data.frame(result_2020$tangency, period = "2020 COVID Crash")
)

period_colors <- c(
  "Full Period (2004–2023)" = "#2C3E7A",
  "2008 Financial Crisis"   = "#C0392B",
  "2020 COVID Crash"        = "#E67E22"
)

p_stress <- ggplot() +

  geom_path(data = all_frontiers,
            aes(x = vol, y = ret, color = period),
            linewidth = 1.1) +

  geom_point(data = tangency_points,
             aes(x = vol, y = ret, color = period),
             size = 4, shape = 18) +

  geom_text(data = tangency_points,
            aes(x = vol, y = ret, label = paste("Tangency\n", period),
                color = period),
            hjust = -0.08, size = 2.8, fontface = "bold", show.legend = FALSE) +

  scale_color_manual(values = period_colors, name = "Period") +

  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +

  labs(
    title    = "Efficient Frontier Shift Under Market Stress",
    subtitle = "Comparing Full Period vs. 2008 Financial Crisis vs. 2020 COVID Crash",
    x        = "Annualized Volatility (σ)",
    y        = "Annualized Expected Return (μ)",
    caption  = "Data: Yahoo Finance via quantmod. Optimization: quadprog (Markowitz QP)."
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(color = "grey40", size = 11),
    plot.caption     = element_text(color = "grey50", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

ggsave("outputs/plots/stress_test.png", p_stress,
       width = 11, height = 7, dpi = 300)

cat("Saved: outputs/plots/stress_test.png\n")

# ── 6. Save stress test results ───────────────────────────────────────────────

saveRDS(list(full = result_full,
             crisis_2008 = result_2008,
             crisis_2020 = result_2020),
        "outputs/data/stress_results.rds")

cat("Stress test complete.\n")