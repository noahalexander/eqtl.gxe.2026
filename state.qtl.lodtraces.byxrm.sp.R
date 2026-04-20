library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

# ---------------- helpers ----------------

prep_cc_long <- function(path_cc_lod, path_state_pheno_lods, time_label) {
  cc <- readRDS(path_cc_lod)
  cc <- as.data.frame(cc)

  cc$state_idx <- seq_len(nrow(cc))

  cc_long <- cc %>%
    mutate(state = state_idx) %>%
    select(state, everything(), -state_idx) %>%
    pivot_longer(
      cols = -state,
      names_to = "pos_chr",
      values_to = "LOD"
    )

  suppressWarnings({
    pos_num <- as.numeric(cc_long$pos_chr)
  })

  if (all(is.na(pos_num))) {
    cc_long <- cc_long %>%
      group_by(state) %>%
      mutate(pos = row_number()) %>%
      ungroup()
  } else {
    cc_long$pos <- pos_num
  }

  roman_states <- c("I","II","III","IV","V","VI","VII","VIII","IX","X")

  cc_long <- cc_long %>%
    mutate(
      state_label = ifelse(state <= length(roman_states),
                           roman_states[state],
                           paste0("S", state)),
      series = paste(time_label, state_label)
    )

  st <- readRDS(path_state_pheno_lods)$LOD
  st <- as.data.frame(t(st))
  st$marker <- rownames(st)
  st$pos <- seq_len(nrow(st))

  chr_map <- st %>%
    transmute(pos,
              chr = str_extract(marker, "^[^_]+")) %>%
    arrange(pos)

  rle_chr <- rle(chr_map$chr)
  cuts <- c(0, cumsum(rle_chr$lengths))
  breaks <- cuts[-length(cuts)]
  labels <- rle_chr$values

  list(data = cc_long,
       breaks = breaks,
       chr_labels = labels)
}

pub_theme <- function() {
  theme_classic(base_size = 18) +
    theme(
      axis.title.y = element_text(size = 22),
      axis.text.x  = element_text(angle = 90,
                                  vjust = 0.5,
                                  hjust = 1,
                                  size = 14),
      axis.text.y  = element_text(size = 16),

      panel.border = element_rect(color = "black",
                                  fill = NA,
                                  size = 1),

      legend.position = "right",
      legend.background =
        element_rect(color = "black",
                     fill = "white",
                     size = 0.8),
      legend.box.background =
        element_rect(color = "black",
                     fill = NA,
                     size = 0.8),

      legend.key = element_blank(),
      legend.text = element_text(size = 16),

      plot.margin = margin(15, 25, 15, 15)
    )
}

# ---------------- t0 ----------------

t0 <- prep_cc_long(
  "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/cell_cycle_assignment_LOD.RDS",
  "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/state_pheno_LODs.RDS",
  "t0"
)

t0_colors <- c(
  "t0 I"   = "lightgreen",
  "t0 II"  = "darkgreen",
  "t0 III" = "lightblue",
  "t0 IV"  = "darkblue",
  "t0 V"   = "black",
  "t0 VI"  = "darkred"
)

p_t0 <- ggplot(t0$data,
               aes(x = pos,
                   y = LOD,
                   color = series)) +
  geom_line(size = 0.7) +
  geom_hline(yintercept = 4,
             size = 0.6,
             color = "grey40") +
  geom_vline(xintercept = t0$breaks,
             linetype = "dashed",
             size = 0.4,
             color = "grey75") +
  scale_color_manual(values = t0_colors,
                     name = NULL) +
  scale_x_continuous(breaks = t0$breaks,
                     labels = t0$chr_labels,
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = NULL,
       y = "LOD") +
  pub_theme()

# ---------------- t10 ----------------

t10 <- prep_cc_long(
  "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/cell_cycle_assignment_LOD.RDS",
  "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/state_pheno_LODs.RDS",
  "t10"
)

t10_colors <- c(
  "t10 I"   = "lightgreen",
  "t10 II"  = "darkgreen",
  "t10 III" = "lightblue",
  "t10 IV"  = "darkblue",
  "t10 V"   = "black",
  "t10 VI"  = "darkred",
  "t10 VII" = "purple"
)

p_t10 <- ggplot(t10$data,
                aes(x = pos,
                    y = LOD,
                    color = series)) +
  geom_line(size = 0.7) +
  geom_hline(yintercept = 4,
             size = 0.6,
             color = "grey40") +
  geom_vline(xintercept = t10$breaks,
             linetype = "dashed",
             size = 0.4,
             color = "grey75") +
  scale_color_manual(values = t10_colors,
                     name = NULL) +
  scale_x_continuous(breaks = t10$breaks,
                     labels = t10$chr_labels,
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = NULL,
       y = "LOD") +
  pub_theme()

# ---------------- combine with panel labels ----------------

combined_plot <- (p_t0 / p_t10) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"),
        plot.tag.position = c(0.02, 0.98))

combined_plot
