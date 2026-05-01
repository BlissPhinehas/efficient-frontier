# 04_plots.R
# Purpose: Plot the efficient frontier, CML, and individual asset points.

library(ggplot2)
library(zoo)

# ── 1. Load results ───────────────────────────────────────────────────────────

efficient_frontier <- readRDS("outputs/data/efficient_frontier.rds")
frontier           <- readRDS("outputs/data/frontier_full.rds")
cml                <- readRDS("outputs/data/cml.rds")
tangency           <- readRDS("outputs/data/tangency.rds")
asset_points       <- readRDS("outputs/data/asset_points.rds")
params             <- readRDS("outputs/data/params.rds")

rf <- params$rf

# ── 2. Main efficient frontier plot ───────────────────────────────────────────

p_frontier <- ggplot() +

  # Full frontier (including inefficient lower half) — faint
  geom_path(data = frontier,
            aes(x = vol, y = ret),
            color = "grey70", linewidth = 0.5, linetype = "dashed") +

  # Efficient frontier — bold
  geom_path(data = efficient_frontier,
            aes(x = vol, y = ret),
            color = "#2C3E7A", linewidth = 1.2) +

  # Capital Market Line
  geom_line(data = cml,
            aes(x = vol, y = ret),
            color = "#C0392B", linewidth = 0.9, linetype = "solid") +

  # Individual asset points
  geom_point(data = asset_points,
             aes(x = vol, y = ret),
             color = "#E67E22", size = 3, shape = 17) +

  # Asset labels
  geom_text(data = asset_points,
            aes(x = vol, y = ret, label = ticker),
            hjust = -0.15, vjust = 0.4, size = 3.2,
            color = "grey30", fontface = "bold") +

  # Tangency portfolio
  geom_point(data = tangency,
             aes(x = vol, y = ret),
             color = "#C0392B", size = 4, shape = 18) +
  annotate("text",
           x = tangency$vol + 0.008, y = tangency$ret,
           label = "Tangency", size = 3.2, color = "#C0392B",
           fontface = "bold", hjust = 0) +

  # Risk-free rate point
  annotate("point", x = 0, y = rf, color = "#C0392B", size = 3) +
  annotate("text",  x = 0.005, y = rf, label = paste0("Rf = ", rf * 100, "%"),
           size = 3, color = "grey30", hjust = 0) +

  # Formatting
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, max(asset_points$vol) * 1.15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Efficient Frontier with Capital Market Line",
    subtitle = "10-Asset Portfolio | 2004–2023 | Long-Only Constraints",
    x        = "Annualized Volatility (σ)",
    y        = "Annualized Expected Return (μ)",
    caption  = "Data: Yahoo Finance via quantmod. Optimization: quadprog (Markowitz QP)."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey40", size = 11),
    plot.caption  = element_text(color = "grey50", size = 9),
    panel.grid.minor = element_blank()
  )

ggsave("outputs/plots/efficient_frontier.png", p_frontier,
       width = 10, height = 7, dpi = 300)

cat("Saved: outputs/plots/efficient_frontier.png\n")

# ── 3. Correlation heatmap ────────────────────────────────────────────────────

log_returns <- readRDS("outputs/data/returns_full.rds")
cor_matrix  <- cor(log_returns)

cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Asset1", "Asset2", "Correlation")

p_corr <- ggplot(cor_df, aes(x = Asset1, y = Asset2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Correlation, 2)),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_gradient2(low  = "#2C3E7A", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(
    title    = "Asset Return Correlation Matrix",
    subtitle = "Daily Log Returns | 2004–2023",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave("outputs/plots/correlation_heatmap.png", p_corr,
       width = 8, height = 7, dpi = 300)

cat("Saved: outputs/plots/correlation_heatmap.png\n")
cat("\nAll plots complete.\n")