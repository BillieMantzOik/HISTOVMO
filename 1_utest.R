# Set working directory and load libraries
setwd("/Volumes/macpor/HISTO_VMO")
library(dplyr)
library(readr)
library(ggplot2)
library(rstatix)
library(ggpubr)

#========
# Load data
#========
histo <- read.csv('./SM2_bfmres.csv')
vars_to_log <- c(
  "collagen_layer_nm","iridophores_area_nm2", "iridophores_thickness_nm", "iridophores_density", 
  "melanophores_area_nm2", "melanophores_thickness_nm", "melanophores_density", 
  "xanthophores_area_nm2", "xanthophore_thickness_nm", "xanthophores_density",
  "xanthtomel", "iritomel", "xanthtoiri", "glands"
)

# morph and species as factors
histo <- histo %>%
  mutate(morph = as.factor(morph),
         species = as.factor(species))
# log1p transform globally
histo <- histo %>%
  mutate(across(all_of(vars_to_log), log1p))
#========
# run Mann-Whitney U test (Wilcoxon rank-sum) per species
#========
run_u_test_by_species <- function(data, species_name, variables) {
  cat("\n============================\n")
  cat("📌 Species:", species_name, "\n")
  cat("============================\n")
  
  df <- data %>%
    filter(species == species_name, morph %in% c("red", "green"))
  %>%
    mutate(across(all_of(variables), log1p))  # log-transform
  
  results <- list()
  
  for (var in variables) {
    cat("\n Variable:", var, "\n")
    
    if (all(is.na(df[[var]]))) {
      cat("Skipped", var, "- all NA\n")
      next
    }
    
    if (length(unique(df$morph[!is.na(df[[var]])])) < 2) {
      cat("Skipped", var, "- less than 2 morph groups\n")
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

species_list <- unique(histo$species)

variables_to_test <- c("collagen_layer_nm","iridophores_area_nm2", "iridophores_thickness_nm", "iridophores_density", 
                       "melanophores_area_nm2", "melanophores_thickness_nm", "melanophores_density", 
                       "xanthophores_area_nm2", "xanthophore_thickness_nm", "xanthophores_density",
                       "xanthtomel", "iritomel", "xanthtoiri", 'glands')
all_utest_results <- list()
for (sp in species_list) {
  res <- run_u_test_by_species(histo, sp, variables_to_test)
  all_utest_results[[sp]] <- res
}

#========
# results
#========
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

#========
# direction of difference (Red > Green or vice versa)
#========
get_median_diff <- function(data, species, variable) {
  df <- data %>%
    filter(species == !!species, morph %in% c("red", "green"))
  
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

species_list <- c("granulifera", "pumilio", "vicentei")
variables_to_test <- c("collagen_layer_nm","iridophores_area_nm2", "iridophores_thickness_nm", "iridophores_density", 
                       "melanophores_area_nm2", "melanophores_thickness_nm", "melanophores_density", 
                       "xanthophores_area_nm2", "xanthophore_thickness_nm", "xanthophores_density",
                       "xanthtomel", "iritomel", "xanthtoiri", "glands")
combinations <- expand.grid(species = species_list, variable = variables_to_test, stringsAsFactors = FALSE)
results_list <- mapply(get_median_diff,
                       species = combinations$species,
                       variable = combinations$variable,
                       MoreArgs = list(data = histo),
                       SIMPLIFY = FALSE)

median_results <- as.data.frame(do.call(rbind, results_list))
median_results <- median_results %>%
  mutate(across(c(red_med, green_med, diff, direction), as.numeric))
