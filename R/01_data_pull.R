# 01_data_pull.R
# Purpose: Pull historical price data, compute log returns,
#          and save clean datasets for optimization and stress testing.

library(quantmod)
library(zoo)
library(dplyr)

# ── 1. Configuration ──────────────────────────────────────────────────────────

tickers <- c("AAPL", "MSFT", "JPM", "XOM", "JNJ",
             "PG", "GLD", "TLT", "SPY", "AMZN")

full_start  <- "2000-01-01"
full_end    <- "2023-12-31"

crisis_2008_start <- "2008-09-01"
crisis_2008_end   <- "2009-03-31"

crisis_2020_start <- "2020-02-01"
crisis_2020_end   <- "2020-04-30"

dir.create("outputs/data",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/plots", recursive = TRUE, showWarnings = FALSE)

# ── 2. Download adjusted closing prices ───────────────────────────────────────

cat("Downloading price data from Yahoo Finance...\n")

price_list <- list()

for (tk in tickers) {
  cat("  Pulling:", tk, "\n")
  tryCatch({
    getSymbols(tk, src = "yahoo", from = full_start, to = full_end,
               auto.assign = TRUE, warnings = FALSE)
    price_list[[tk]] <- Ad(get(tk))  # adjusted close only
  }, error = function(e) {
    cat("  WARNING: Failed to pull", tk, "-", conditionMessage(e), "\n")
  })
}

# ── 3. Align by date intersection ─────────────────────────────────────────────
# GLD (est. 2004) and TLT (est. 2002) don't go back to 2000.
# We take the intersection of all available dates so the
# covariance matrix is always computed on identical observations.

prices_merged <- do.call(merge, price_list)
colnames(prices_merged) <- names(price_list)

# Drop any rows where ANY ticker has NA
prices_merged <- na.omit(prices_merged)

cat("\nDate range after alignment:", 
    as.character(index(prices_merged)[1]), "to",
    as.character(index(prices_merged)[nrow(prices_merged)]), "\n")
cat("Total trading days:", nrow(prices_merged), "\n")
cat("Assets included:   ", ncol(prices_merged), "\n")

# ── 4. Compute daily log returns ──────────────────────────────────────────────

log_returns <- diff(log(prices_merged))
log_returns <- na.omit(log_returns)

cat("\nReturns matrix dimensions:", dim(log_returns), "\n")

# ── 5. Carve out crisis windows ───────────────────────────────────────────────

returns_2008 <- window(log_returns,
                       start = as.Date(crisis_2008_start),
                       end   = as.Date(crisis_2008_end))

returns_2020 <- window(log_returns,
                       start = as.Date(crisis_2020_start),
                       end   = as.Date(crisis_2020_end))

cat("\n2008 crisis window:", nrow(returns_2008), "trading days\n")
cat("2020 crisis window:", nrow(returns_2020), "trading days\n")

# ── 6. Save to disk ───────────────────────────────────────────────────────────

saveRDS(log_returns,    "outputs/data/returns_full.rds")
saveRDS(returns_2008,  "outputs/data/returns_2008.rds")
saveRDS(returns_2020,  "outputs/data/returns_2020.rds")
saveRDS(prices_merged, "outputs/data/prices_full.rds")

cat("\nAll datasets saved to outputs/data/\n")
cat("Data pull complete.\n")