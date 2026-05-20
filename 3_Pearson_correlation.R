# Script by V. Mantzana-Oikonomaki for CCA-analysis (CLEANED)

library(tidyverse)
library(Hmisc)
library(corrplot)
library(ggplot2)

setwd("/Volumes/macpor/HISTO_VMO/")

# ============================================================
# 1. Load & clean data
# ============================================================
histo <- read.csv('./tosubmit/SM2.csv') %>%
  mutate(
    species = trimws(as.character(species)),
    morph   = trimws(as.character(morph))
  ) %>%
  filter(morph %in% c("red", "green"))

# ============================================================
# 2. Compute relative chromatophore investment
# ============================================================
histo <- histo %>%
  mutate(
    total_chrom_area = xanthophores_area_px +
      melanophores_area_px +
      iridophores_area_px,
    
    rel_xanth = ifelse(total_chrom_area > 0, xanthophores_area_px / total_chrom_area, NA),
    rel_melan = ifelse(total_chrom_area > 0, melanophores_area_px / total_chrom_area, NA),
    rel_irid  = ifelse(total_chrom_area > 0, iridophores_area_px / total_chrom_area, NA)
  )

# ============================================================
# 3. Define variable sets (CLEANED)
# ============================================================

spectral_vars <- c("B3", "S9", "S1V", "S1G")

histo_vars <- c(
  "collagen_layer_nm",
  "iridophores_area_nm2",
  "iridophores_thickness_nm",
  "melanophores_area_nm2",
  "melanophores_thickness_nm",
  "xanthophores_area_nm2",
  "xanthophore_thickness_nm",
  "rel_xanth",
  "rel_melan",
  "rel_irid"
)

stopifnot(all(c(spectral_vars, histo_vars) %in% colnames(histo)))

# ============================================================
# 4. Remove rows with NA in either group
# ============================================================
histo_clean <- histo %>% drop_na(all_of(c(spectral_vars, histo_vars)))

# ============================================================
# 5. Log-transform histological variables (avoid log(0))
# ============================================================
histo_clean <- histo_clean %>%
  mutate(across(all_of(histo_vars), ~ log(.x + 1)))

# ============================================================
# 6. Split by morph
# ============================================================
histo_red   <- histo_clean %>% filter(morph == "red")
histo_green <- histo_clean %>% filter(morph == "green")

# ============================================================
# 7. Correlation function
# ============================================================
get_correlation_matrix <- function(data, x_vars, y_vars) {
  x <- data[, x_vars, drop = FALSE]
  y <- data[, y_vars, drop = FALSE]
  keep <- complete.cases(x, y)
  
  if (sum(keep) < 5) stop("Not enough complete cases for correlation analysis.")
  
  cor_matrix <- rcorr(as.matrix(x[keep, ]), as.matrix(y[keep, ]))
  
  list(
    cor = cor_matrix$r[x_vars, y_vars],
    p   = cor_matrix$P[x_vars, y_vars]
  )
}

# ============================================================
# 8. Run correlations
# ============================================================
cor_red   <- get_correlation_matrix(histo_red, spectral_vars, histo_vars)
cor_green <- get_correlation_matrix(histo_green, spectral_vars, histo_vars)

# ============================================================
# 9. Plot correlation matrices
# ============================================================

corrplot(cor_red$cor,
         p.mat = cor_red$p,
         method = "color",
         type = "full",
         sig.level = 0.05,
         insig = "blank",
         tl.col = "black",
         title = "Red Frogs (Log-Transformed Histology)",
         mar = c(0, 0, 2, 0))

corrplot(cor_green$cor,
         p.mat = cor_green$p,
         method = "color",
         type = "full",
         sig.level = 0.05,
         insig = "blank",
         tl.col = "black",
         title = "Green Frogs (Log-Transformed Histology)",
         mar = c(0, 0, 2, 0))

# ============================================================
# 10. Extract significant correlations
# ============================================================
extract_significant_correlations <- function(cor_list, color_label, r_thresh = 0.4, p_thresh = 0.05) {
  
  cor_df <- as.data.frame(as.table(cor_list$cor))
  p_df   <- as.data.frame(as.table(cor_list$p))
  
  df <- merge(cor_df, p_df, by = c("Var1", "Var2"))
  colnames(df) <- c("Spectral", "Histological", "Correlation", "p_value")
  
  df %>%
    filter(abs(Correlation) >= r_thresh, p_value <= p_thresh) %>%
    mutate(Color = color_label)
}

red_df   <- extract_significant_correlations(cor_red, "Red")
green_df <- extract_significant_correlations(cor_green, "Green")

combined_df <- bind_rows(red_df, green_df)

# ============================================================
# 11. Barplot of significant correlations
# ============================================================
ggplot(combined_df,
       aes(x = interaction(Spectral, Histological),
           y = Correlation,
           fill = Color)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Red" = "#BF193A", "Green" = "#57F257")) +
  labs(title = "Significant Spectral–Histological Correlations (Log-Transformed Histology)",
       x = "Spectral ~ Histological Pair",
       y = "Pearson Correlation",
       fill = "Frog Color") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

ggsave('pearson_cor_sphis_log_v2.pdf', width = 10, height = 8, dpi = 600)
