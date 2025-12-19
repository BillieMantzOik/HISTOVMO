#Script by V. Mantzana-Oikonomaki for CCA-analysis
# Load libraries
library(tidyverse)
library(Hmisc)
library(corrplot)
library(ggplot2)
setwd("/Volumes/macpor/HISTO_VMO/")
histo <- read.csv('./tosubmit/SM2_bfmres.csv')
# Define variable sets
spectral_vars <- c("B3", "S9", "S1V", "S1G")
histo_vars <- c("collagen_layer_nm", "iridophores_area_nm2", "iridophores_thickness_nm", 
                "iridophores_density", "melanophores_area_nm2", "melanophores_thickness_nm", 
                "melanophores_density", "xanthophores_area_nm2", "xanthophore_thickness_nm",
                "xanthophores_density", "xanthtomel", "iritomel", "xanthtoiri", "glands")
# Check variable presence
stopifnot(all(c(spectral_vars, histo_vars) %in% colnames(histo)))

# Remove rows with NA in either group
histo_clean <- histo %>% drop_na(all_of(c(spectral_vars, histo_vars)))

# Log-transform histological variables (add +1 to avoid log(0))
histo_clean <- histo_clean %>%
  mutate(across(all_of(histo_vars), ~ log(.x + 1)))

# Split by color morph
histo_red <- histo_clean %>% filter(morph == "red")
histo_green <- histo_clean %>% filter(morph == "green")

#Define correlation function
get_correlation_matrix <- function(data, x_vars, y_vars) {
  x <- data[, x_vars, drop = FALSE]
  y <- data[, y_vars, drop = FALSE]
  keep <- complete.cases(x, y)
  
  if (sum(keep) < 5) stop("Not enough complete cases for correlation analysis.")
  
  cor_matrix <- rcorr(as.matrix(x[keep, ]), as.matrix(y[keep, ]))
  list(cor = cor_matrix$r[x_vars, y_vars], p = cor_matrix$P[x_vars, y_vars])
}
# Run correlations
cor_red <- get_correlation_matrix(histo_red, spectral_vars, histo_vars)
cor_green <- get_correlation_matrix(histo_green, spectral_vars, histo_vars)

# Plot correlation matrices
# Red frogs
corrplot(cor_red$cor,
         p.mat = cor_red$p,
         method = "color",
         type = "full",
         sig.level = 0.05,
         insig = "blank",
         tl.col = "black",
         title = "Red Frogs (Log-Transformed Histology)",
         mar = c(0, 0, 2, 0))

# Green frogs
corrplot(cor_green$cor,
         p.mat = cor_green$p,
         method = "color",
         type = "full",
         sig.level = 0.05,
         insig = "blank",
         tl.col = "black",
         title = "Green Frogs (Log-Transformed Histology)",
         mar = c(0, 0, 2, 0))
# Extract and plot significant correlations
extract_significant_correlations <- function(cor_list, color_label, r_thresh = 0.4, p_thresh = 0.05) {
  cor_df <- as.data.frame(as.table(cor_list$cor))
  p_df   <- as.data.frame(as.table(cor_list$p))
  df <- merge(cor_df, p_df, by = c("Var1", "Var2"))
  colnames(df) <- c("Spectral", "Histological", "Correlation", "p_value")
  
  df %>%
    filter(abs(Correlation) >= r_thresh, p_value <= p_thresh) %>%
    mutate(Color = color_label)
}
# Get significant results
red_df   <- extract_significant_correlations(cor_red, "Red")
green_df <- extract_significant_correlations(cor_green, "Green")
combined_df <- bind_rows(red_df, green_df)
