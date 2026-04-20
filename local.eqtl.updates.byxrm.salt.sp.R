# ============================================================
# Full script: build two-panel figure from start to finish
# - Panel A: NaCl t0 vs t30
# - Panel B: SP t0 vs t10
# - One shared legend on the right
# - Legend has no box
# - Legend title is bold, size 15
# - Panel labels A and B in upper left
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ----------------------------
# Read input files
# ----------------------------
df.nt0 <- readRDS("repeat_fine_mapping/combined/A/nacl/t0/cis_only_test_CombinedResults.RDS")
df.nt0 <- df.nt0$combined

df.3 <- readRDS("repeat_fine_mapping/combined/A/nacl/t30/cis_only_test_CombinedResults.RDS")
df.3 <- df.3$combined

df.spt0 <- readRDS("repeat_fine_mapping/combined/A/sp/t0/cis_only_test_CombinedResults.RDS")
df.spt0 <- df.spt0$combined

df.t10 <- readRDS("repeat_fine_mapping/combined/A/sp/t10/cis_only_test_CombinedResults.RDS")
df.t10 <- df.t10$combined

# ----------------------------
# Build merged data for panel A
# NaCl: t0 vs t30
# ----------------------------
dl_nt <- merge(df.nt0, df.3, by = "transcript", all = TRUE)

dl_nt$FDR.x[is.na(dl_nt$FDR.x)] <- 1
dl_nt$FDR.y[is.na(dl_nt$FDR.y)] <- 1

dl_nt$vec <- case_when(
  dl_nt$FDR.x <= 0.05 & dl_nt$FDR.y <= 0.05 ~ "Both timepoints",
  dl_nt$FDR.x >  0.05 & dl_nt$FDR.y <= 0.05 ~ "Perturbation only",
  dl_nt$FDR.x <= 0.05 & dl_nt$FDR.y >  0.05 ~ "Baseline only",
  TRUE                                      ~ "Neither"
)

# ----------------------------
# Build merged data for panel B
# SP: t0 vs t10
# ----------------------------
dl_sp <- merge(df.spt0, df.t10, by = "transcript", all = TRUE)

dl_sp$FDR.x[is.na(dl_sp$FDR.x)] <- 1
dl_sp$FDR.y[is.na(dl_sp$FDR.y)] <- 1

dl_sp$vec <- case_when(
  dl_sp$FDR.x <= 0.05 & dl_sp$FDR.y <= 0.05 ~ "Both timepoints",
  dl_sp$FDR.x >  0.05 & dl_sp$FDR.y <= 0.05 ~ "Perturbation only",
  dl_sp$FDR.x <= 0.05 & dl_sp$FDR.y >  0.05 ~ "Baseline only",
  TRUE                                      ~ "Neither"
)

# ----------------------------
# Force identical factor levels
# so the shared legend behaves correctly
# ----------------------------
common_levels <- c("Both timepoints", "Perturbation only", "Baseline only", "Neither")

dl_nt$vec <- factor(dl_nt$vec, levels = common_levels)
dl_sp$vec <- factor(dl_sp$vec, levels = common_levels)

# ----------------------------
# Shared color scale
# ----------------------------
sig_scale <- scale_color_manual(
  name   = "Significance\n(FDR \u2264 0.05)",
  breaks = common_levels,
  limits = common_levels,
  drop   = FALSE,
  values = c(
    "Both timepoints"   = "#D55E00",
    "Perturbation only" = "#0072B2",
    "Baseline only"     = "#009E73",
    "Neither"           = "grey70"
  )
)

# ----------------------------
# Panel A
# ----------------------------
pA <- ggplot(dl_nt, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point(size = 1.2, alpha = 0.6) +
  sig_scale +
  labs(
    x = "log(LOD) at t0 time point",
    y = "log(LOD) at t30 time point"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "A",
    hjust = -0.5, vjust = 1.5,
    size = 6, fontface = "bold"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

# ----------------------------
# Panel B
# ----------------------------
pB <- ggplot(dl_sp, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point(size = 1.2, alpha = 0.6) +
  sig_scale +
  labs(
    x = "log(LOD) at t0 time point",
    y = "log(LOD) at t10 time point"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "B",
    hjust = -0.5, vjust = 1.5,
    size = 6, fontface = "bold"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 11, color = "black"),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# ----------------------------
# Combine panels with one shared legend
# ----------------------------
final_plot <- (pA + pB) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12)
  )

# Show plot
final_plot

# ----------------------------
# Optional export
# Uncomment if needed
# ----------------------------
# ggsave(
#   filename = "two_panel_LOD_plot.pdf",
#   plot = final_plot,
#   width = 12,
#   height = 6,
#   units = "in"
# )
