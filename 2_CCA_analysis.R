# Script by V. Mantzana-Oikonomaki for CCA-analysis
library(vegan)      
library(tidyverse)   
library(ggfortify)   
library(ggrepel)     
library(ggplot2)

setwd("/Volumes/macpor/HISTO_VMO/")
histo <- read.csv('./tosubmit/SM2.csv')
histo <- histo %>% filter(morph %in% c("red", "green"))

# =====================================================================
# PCA Analysis ---------------------------------------------------------
# =====================================================================

spectral_cols <- c(
  "B1", "B2", "B3",
  "S1U", "S1V", "S1B", "S1G", "S1Y", "S1R",
  "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
  "H1", "H2", "H3", "H4", "H5"
)

frog_spectral_data <- histo[, c("id", "morph", "species", spectral_cols)]
frog_spectral_data_clean <- frog_spectral_data[complete.cases(frog_spectral_data[, spectral_cols]), ]

# PCA computation
pca_spec <- prcomp(frog_spectral_data_clean[, spectral_cols], scale. = TRUE)

# Scores and loadings
scores <- as.data.frame(pca_spec$x)
scores$id <- frog_spectral_data_clean$id
scores$morph <- frog_spectral_data_clean$morph
scores$species <- frog_spectral_data_clean$species

loadings <- as.data.frame(pca_spec$rotation[, 1:2])
loadings$variable <- rownames(loadings)

# Scale arrows by loading strength
scale_factor <- 5
loadings$PC1 <- loadings$PC1 * scale_factor
loadings$PC2 <- loadings$PC2 * scale_factor

# PCA biplot
ggplot() +
  geom_point(data = scores, aes(x = PC1, y = PC2, color = morph, shape = species),
             size = 3, alpha = 0.8) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "gray40", linewidth = 1) +
  geom_text_repel(data = loadings, aes(x = PC1, y = PC2, label = variable),
                  color = "black", size = 5) +
  scale_color_manual(values = c("#57F257", "#BF193A")) +
  scale_shape_manual(values = c(17, 16, 15)) +
  theme_minimal(base_size = 20) +
  labs(x = "PC1", y = "PC2", color = "Morph", shape = "Species",
       title = "PCA Biplot: Samples and Loadings") +
  theme(legend.position = "right")

ggsave("pca_spectral_V2.pdf", width = 8, height = 6, dpi = 900)

# =====================================================================
# CCA Analysis ---------------------------------------------------------
# =====================================================================

# Response variables (spectral)
response_vars <- c("B3", "S9", "S1V", "S1G", "S1R")

# *** Corrected histology variables ***
explanatory_vars <- c(
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

Y <- histo %>% select(all_of(response_vars))
X <- histo %>% select(all_of(explanatory_vars))

# Remove rows with missing values
common_rows <- complete.cases(Y, X)
Y <- Y[common_rows, ]
X <- X[common_rows, ]
histo_sub <- histo[common_rows, ]

# Log-transform explanatory variables ONLY (correct for CCA)
X_log <- X %>% mutate(across(everything(), ~ log(.x + 1)))

# Run CCA
cca_model <- cca(Y ~ ., data = X_log)

# Summary and significance tests
summary(cca_model)
anova(cca_model, by = "term")
anova(cca_model, by = "axis")

# Extract scores for plotting
sites <- as.data.frame(scores(cca_model, display = "sites", scaling = 2))
sites$id <- histo_sub$id
sites$morph <- histo_sub$morph
sites$species <- histo_sub$species

variables <- as.data.frame(scores(cca_model, display = "bp", scaling = 2))
variables$label <- rownames(variables)

# Strong variables (same threshold as your script)
variables_strong <- subset(variables, abs(CCA1) > 0.3 | abs(CCA2) > 0.3)

# CCA biplot
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2, color = morph, shape = species),
             size = 3, alpha = 0.8) +
  geom_segment(data = variables,
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30", linewidth = 0.6) +
  geom_text_repel(data = variables_strong,
                  aes(x = CCA1, y = CCA2, label = label),
                  size = 6, max.overlaps = 100) +
  scale_color_manual(values = c("#57F257", "#BF193A")) +
  scale_shape_manual(values = c(17, 16, 15)) +
  theme_minimal(base_size = 20)

ggsave('cca_sphisto_V2.pdf', width = 8, height = 6, dpi = 900)

# =====================================================================
# Correlation of log-transformed variables with CCA1 -------------------
# =====================================================================

cca_scores <- scores(cca_model, display = "sites", scaling = 2)
CCA1 <- cca_scores[, "CCA1"]

cor_results <- sapply(X_log, function(v) cor(v, CCA1, method = "pearson"))

cor_df <- data.frame(
  Variable = names(cor_results),
  Correlation = cor_results
) %>%
  arrange(desc(Correlation))

ggplot(cor_df, aes(x = reorder(Variable, Correlation), y = Correlation)) +
  geom_col(fill = "darkcyan") +
  geom_hline(yintercept = 0, linewidth = 1.2) +
  
  # ---- r-values inside bars ----
geom_text(
  aes(label = sprintf("%.2f", Correlation)),
  color = "white",
  size = 5,
  hjust = ifelse(cor_df$Correlation > 0, 1.1, -0.1)  # inside bar on correct side
) +
  
  coord_flip() +
  theme_minimal(base_size = 20) +
  labs(
    title = "Correlation of Variables with CCA1",
    x = "Histological Variable",
    y = "Pearson Correlation"
  )

ggsave("Figure6C_CCA1_loadings_V2.pdf", width = 8, height = 6, dpi = 900)


# =====================================================================
# Export loadings ------------------------------------------------------
# =====================================================================

resp_loadings_s1 <- as.data.frame(scores(cca_model, display = "species", scaling = 1))
resp_loadings_s2 <- as.data.frame(scores(cca_model, display = "species", scaling = 2))

resp_loadings_s1$Variable <- rownames(resp_loadings_s1)
resp_loadings_s2$Variable <- rownames(resp_loadings_s2)

write.csv(resp_loadings_s1, "cca_response_loadings_scaling1_V2.csv", row.names = FALSE)
write.csv(resp_loadings_s2, "cca_response_loadings_scaling2_V2.csv", row.names = FALSE)

env_loadings_s1 <- as.data.frame(scores(cca_model, display = "bp", scaling = 1))
env_loadings_s2 <- as.data.frame(scores(cca_model, display = "bp", scaling = 2))

env_loadings_s1$Variable <- rownames(env_loadings_s1)
env_loadings_s2$Variable <- rownames(env_loadings_s2)

write.csv(env_loadings_s1, "cca_env_loadings_scaling1_V2.csv", row.names = FALSE)
write.csv(env_loadings_s2, "cca_env_loadings_scaling2_V2.csv", row.names = FALSE)

canon_coeff <- as.data.frame(coef(cca_model))
canon_coeff$Variable <- rownames(canon_coeff)

write.csv(canon_coeff, "cca_canonical_coefficients_V2.csv", row.names = FALSE)
