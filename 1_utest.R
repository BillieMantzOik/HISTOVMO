# ============================================
# 0. Setup
# ============================================
setwd("/Users/vasiliki/desktop/HISTO_VMO/tosubmit")

library(dplyr)
library(readr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(janitor)

# ============================================
# 1. Load data + CLEAN species/morph + compute relative investment
# ============================================
histo <- read.csv('./SM2.csv') %>%
  janitor::clean_names() %>%   # clean column names
  mutate(
    species = trimws(as.character(species)),   # CLEAN species names
    morph   = trimws(as.character(morph))      # CLEAN morph names
  )

# Compute total chromatophore area + relative investment
histo <- histo %>%
  mutate(
    total_chrom_area = xanthophores_area_px +
      melanophores_area_px +
      iridophores_area_px,
    
    rel_xanth = ifelse(total_chrom_area > 0, xanthophores_area_px / total_chrom_area, NA),
    rel_melan = ifelse(total_chrom_area > 0, melanophores_area_px / total_chrom_area, NA),
    rel_irid  = ifelse(total_chrom_area > 0, iridophores_area_px / total_chrom_area, NA)
  )

# Create species list AFTER cleaning
species_list <- unique(histo$species)

# ============================================
# 2. Mann–Whitney U test function
# ============================================
run_u_test_by_species <- function(data, species_name, variables) {
  cat("\n============================\n")
  cat("📌 Species:", species_name, "\n")
  cat("============================\n")
  
  df <- data %>%
    filter(species == species_name, morph %in% c("red", "green"))
  
  results <- list()
  
  for (var in variables) {
    cat("\n🔹 Variable:", var, "\n")
    
    if (all(is.na(df[[var]]))) {
      cat("⚠️ Skipped", var, "- all NA\n")
      next
    }
    
    if (length(unique(df$morph[!is.na(df[[var]])])) < 2) {
      cat("⚠️ Skipped", var, "- less than 2 morph groups\n")
      next
    }
    
    test_result <- wilcox.test(
      x = df[[var]][df$morph == "red"],
      y = df[[var]][df$morph == "green"],
      exact = FALSE,
      paired = FALSE
    )
    
    results[[var]] <- test_result
  }
  
  return(results)
}

# ============================================
# 3. Variables to test
# ============================================
variables_to_test <- c(
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

# Run U-tests
all_utest_results <- list()
for (sp in species_list) {
  res <- run_u_test_by_species(histo, sp, variables_to_test)
  all_utest_results[[sp]] <- res
}

# ============================================
# 4. Extract U-test results
# ============================================
extract_u_test_results <- function(all_utest_results) {
  df_list <- list()
  
  for (species_name in names(all_utest_results)) {
    species_results <- all_utest_results[[species_name]]
    
    for (var_name in names(species_results)) {
      test <- species_results[[var_name]]
      
      df_list[[paste(species_name, var_name, sep = "_")]] <- data.frame(
        species = species_name,
        variable = var_name,
        p.value = test$p.value,
        statistic = test$statistic
      )
    }
  }
  
  results_df <- do.call(rbind, df_list)
  rownames(results_df) <- NULL
  
  results_df$significant <- ifelse(results_df$p.value < 0.05, "Significant", "Not Significant")
  
  return(results_df)
}

u_test_df <- extract_u_test_results(all_utest_results)

# ============================================
# 5. Compute direction of difference
# ============================================
get_median_diff <- function(data, species, variable) {
  df <- data %>%
    filter(species == species, morph %in% c("red", "green"))
  
  red_med <- median(df[[variable]][df$morph == "red"], na.rm = TRUE)
  green_med <- median(df[[variable]][df$morph == "green"], na.rm = TRUE)
  
  return(c(
    species = species,
    variable = variable,
    red_med = red_med,
    green_med = green_med,
    diff = red_med - green_med,
    direction = sign(red_med - green_med)
  ))
}

combinations <- expand.grid(
  species = species_list,
  variable = variables_to_test,
  stringsAsFactors = FALSE
)

results_list <- mapply(
  get_median_diff,
  species = combinations$species,
  variable = combinations$variable,
  MoreArgs = list(data = histo),
  SIMPLIFY = FALSE
)

median_results <- as.data.frame(do.call(rbind, results_list)) %>%
  mutate(across(c(red_med, green_med, diff, direction), as.numeric))

# ============================================
# 6. Merge U-test results with medians
# ============================================
u_test_df <- u_test_df %>%
  left_join(median_results, by = c("species", "variable")) %>%
  mutate(
    direction_label = ifelse(direction == 1, "Red > Green",
                             ifelse(direction == -1, "Green > Red", "Equal"))
  )

# ============================================
# 7. Heatmap of significant differences
# ============================================
heatmap_df <- u_test_df %>% filter(significant == "Significant")

ggplot(heatmap_df, aes(x = variable, y = species, fill = direction_label)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Red > Green" = "#BF193A",
                               "Green > Red" = "#57F257",
                               "Equal" = "grey80")) +
  labs(title = "Direction of Significant Differences (Red vs. Green)",
       x = "Histological Variable",
       y = "Species",
       fill = "Direction") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("utest_bfm_last.pdf", dpi = 600, width = 10, height = 7)

# ============================================
# 8. Summary tables
# ============================================
get_summary_stats <- function(data, species, variable) {
  df <- data %>%
    filter(species == species, morph %in% c("red", "green"))
  
  df %>%
    group_by(species, morph) %>%
    summarise(
      median = median(.data[[variable]], na.rm = TRUE),
      q1 = quantile(.data[[variable]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[variable]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

build_summary_table <- function(data, u_test_df) {
  summary_rows <- list()
  sig_df <- u_test_df %>% filter(significant == "Significant")
  
  for (i in 1:nrow(sig_df)) {
    species <- sig_df$species[i]
    variable <- sig_df$variable[i]
    
    stat_row <- get_summary_stats(data, species, variable)
    red_vals <- stat_row %>% filter(morph == "red")
    green_vals <- stat_row %>% filter(morph == "green")
    
    summary_rows[[i]] <- data.frame(
      species = species,
      variable = variable,
      red_median_iqr = sprintf("%.2f (%.2f–%.2f)", red_vals$median, red_vals$q1, red_vals$q3),
      green_median_iqr = sprintf("%.2f (%.2f–%.2f)", green_vals$median, green_vals$q1, green_vals$q3),
      U_statistic = sig_df$statistic[i],
      p_value = sig_df$p.value[i],
      direction = sig_df$direction_label[i]
    )
  }
  
  return(do.call(rbind, summary_rows))
}

summary_table <- build_summary_table(histo, u_test_df)
write.csv(summary_table, "significant_results_table_bfm_nolog_v2.csv", row.names = FALSE)

# Full summary table
build_full_summary_table <- function(data, u_test_df) {
  summary_rows <- list()
  
  for (i in 1:nrow(u_test_df)) {
    species <- u_test_df$species[i]
    variable <- u_test_df$variable[i]
    
    stat_row <- get_summary_stats(data, species, variable)
    red_vals <- stat_row %>% filter(morph == "red")
    green_vals <- stat_row %>% filter(morph == "green")
    
    summary_rows[[i]] <- data.frame(
      species = species,
      variable = variable,
      red_median_iqr = ifelse(nrow(red_vals) == 0, NA,
                              sprintf("%.2f (%.2f–%.2f)", red_vals$median, red_vals$q1, red_vals$q3)),
      green_median_iqr = ifelse(nrow(green_vals) == 0, NA,
                                sprintf("%.2f (%.2f–%.2f)", green_vals$median, green_vals$q1, green_vals$q3)),
      U_statistic = u_test_df$statistic[i],
      p_value = u_test_df$p.value[i],
      significant = u_test_df$significant[i],
      direction = u_test_df$direction_label[i]
    )
  }
  
  return(do.call(rbind, summary_rows))
}

full_summary_table <- build_full_summary_table(histo, u_test_df)
write.csv(full_summary_table, "full_results_table_bfm_nolog_v2.csv", row.names = FALSE)
