### Install packages
install.packages(c("tidyverse","ggplot2","patchwork","lme4","lmerTest","emmeans","performance","car", "ggcorrplot"))
 
library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)
library(ggcorrplot)

### LOAD DATA

setwd("C:/Users/RebeccaPedler/OneDrive - Yumbah/Documents/R&D/Industry PhD/Trials/Commercial trial/R_datasets")
df_raw <- read.csv("individual_abalone_data.csv")
str(df_raw)

# Make tank and diet factors
df_raw <- df_raw |>
  mutate(
    tank = factor(tank),
    diet = factor(diet, levels = c("control", "ulva", "wakame")),
 
    log_weight = log(weight_g)
  )

## Check for missing data
missing_summary <- df_raw |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  mutate(pct_missing = round(n_missing / nrow(df) * 100, 2)) |>
  filter(n_missing > 0)

print(missing_summary) # dataframe is empty

# Filter out biologically implausible values (abalone below 10g or above 150g) 
n_before <- nrow(df_raw)

df <- df_raw |>
  filter(
    length_mm > 30,
    weight_g >= 10,
    weight_g  <= 150
  )

n_removed <- n_before - nrow(df)
cat(sprintf("  Removed %d rows (%d remaining)\n", n_removed, nrow(df)))

# Create colour palette for plots
diet_cols <- c(
  "control" = "#8DB4C8",   
  "ulva"    = "#6DAA6E",   
  "wakame"  = "#C8844A"    
)

### CLEAN DATA

## Check sample sizes
df |>
  count(diet, name = "n_individuals") |>
  as.data.frame() |>
  print(row.names = FALSE)
 
# Observations per tank × diet
df |>
  count(tank, diet, name = "n") |>
  arrange(diet, tank) |>
  as.data.frame() |>
  print(row.names = FALSE)
 
## Detect outliers, clean and generate numerical summarise

# Detect outliers in length_mm only
flag_outliers <- function(data, var) {
  data |>
    group_by(tank) |>
    mutate(
      tank_mean = mean(.data[[var]], na.rm = TRUE),
      tank_sd   = sd(.data[[var]],   na.rm = TRUE),
      z_within  = (.data[[var]] - tank_mean) / tank_sd,
      outlier   = abs(z_within) > 3
    ) |>
    ungroup() |>
    filter(outlier) |>
    mutate(
      direction = if_else(z_within < 0, "small", "large")
    ) |>
    dplyr::select(tank, diet, all_of(var), z_within, direction) |>
    mutate(across(where(is.numeric), ~ round(.x, 3))) |>
    arrange(z_within)
}

# Inspect for length_mm
out <- flag_outliers(df, "length_mm")
cat(sprintf("\n  Within-tank outliers (|z|>3) for length_mm: %d rows\n", nrow(out)))if (nrow(out) > 0) {
  small <- filter(out, direction == "small")
  large <- filter(out, direction == "large")
  
  if (nrow(small) > 0) {
    cat(sprintf("    -- Small outliers (z < -3): %d rows\n", nrow(small)))
    print(as.data.frame(small), row.names = FALSE)
  }
  if (nrow(large) > 0) {
    cat(sprintf("    -- Large outliers (z >  3): %d rows\n", nrow(large)))
    print(as.data.frame(large), row.names = FALSE)
  }
}

# Inspect for weight_g
out <- flag_outliers(df, "weight_g")
cat(sprintf("\n  Within-tank outliers (|z|>3) for weight_g: %d rows\n", nrow(out)))if (nrow(out) > 0) {
  small <- filter(out, direction == "small")
  large <- filter(out, direction == "large")
  
  if (nrow(small) > 0) {
    cat(sprintf("    -- Small outliers (z < -3): %d rows\n", nrow(small)))
    print(as.data.frame(small), row.names = FALSE)
  }
  if (nrow(large) > 0) {
    cat(sprintf("    -- Large outliers (z >  3): %d rows\n", nrow(large)))
    print(as.data.frame(large), row.names = FALSE)
  }
}

### INSPECT DATA

## Histograms

# Helper: histogram + density coloured by diet
hist_diet <- function(var, xlab, binwidth = NULL) {
  p <- ggplot(df, aes(x = .data[[var]], fill = diet, colour = diet)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins     = 45,
      alpha    = 0.45,
      position = "identity",
      linewidth = 0.2
    ) +
    geom_density(alpha = 0, linewidth = 0.7) +
    scale_fill_manual(values = diet_cols,   name = "Diet") +
    scale_colour_manual(values = diet_cols, name = "Diet") +
    labs(x = xlab, y = "Density") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
  p
}
 
# Helper: boxplot by tank, coloured by diet
box_tank <- function(var, ylab) {
  ggplot(df, aes(x = tank, y = .data[[var]], fill = diet)) +
    geom_boxplot(
      outlier.size  = 0.7,
      outlier.alpha = 0.4,
      colour        = "grey30",
      linewidth     = 0.35,
      alpha         = 0.75
    ) +
    scale_fill_manual(values = diet_cols, name = "Diet") +
    labs(x = "Tank", y = ylab) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}
 
# length
p_len_hist <- hist_diet("length_mm", "Final shell length (mm)") + ggtitle("Shell length distribution by diet")
p_len_hist
 
p_len_box  <- box_tank("length_mm", "Final shell length (mm)") + ggtitle("Shell length by tank (coloured by diet)")
p_len_box

# weight (raw and log)
p_wt_hist <- hist_diet("weight_g", "Final weight (g)") + ggtitle("Weight distribution by diet")
p_wt_hist
 
p_wt_log_hist <- hist_diet("log_weight", "log(weight) (g)") + ggtitle("log(Weight) distribution by diet")
p_wt_log_hist
 
p_wt_box  <- box_tank("weight_g", "Final weight (g)") + ggtitle("Weight by tank (coloured by diet)")
p_wt_box

## Numerical summaries

# Define variables
response_vars <- c("length_mm", "weight_g", "condition")

# Create table
for (v in response_vars) {
  cat(sprintf("\n--- %s ---\n", v))
  out <- df |>
    group_by(diet) |>
    summarise(
      n      = n(),
      mean   = round(mean(.data[[v]], na.rm = TRUE), 3),
      sd     = round(sd(.data[[v]],   na.rm = TRUE), 3),
      min    = round(min(.data[[v]],  na.rm = TRUE), 3),
      q25    = round(quantile(.data[[v]], 0.25, na.rm = TRUE), 3),
      median = round(median(.data[[v]], na.rm = TRUE), 3),
      q75    = round(quantile(.data[[v]], 0.75, na.rm = TRUE), 3),
      max    = round(max(.data[[v]],  na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    as.data.frame()
  print(out, row.names = FALSE)
}

## CHECK ALL STARTING CONDITIONS

# Collapse to one row per tank — start values are constant within tank
tank_df <- df |>
  group_by(tank, diet) |>
  summarise(
    start_ABL     = first(start_ABL),
    start_ABW     = first(start_ABW),
    start_density = first(start_density),
    start_biomass = first(start_biomass),
    .groups = "drop"
  )
 
# Print tank-level table
cat("\n  Starting values per tank:\n")
tank_df |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)
 
# Summary by diet
cat("\n  Starting values summarised by diet:\n")
start_vars <- c("start_ABL", "start_ABW", "start_density", "start_biomass")
for (v in start_vars) {
  cat(sprintf("\n--- %s ---\n", v))
  tank_df |>
    group_by(diet) |>
    summarise(
      n      = n(),
      mean   = round(mean(.data[[v]]),   3),
      sd     = round(sd(.data[[v]]),     3),
      min    = round(min(.data[[v]]),    3),
      max    = round(max(.data[[v]]),    3),
      range  = round(max(.data[[v]]) - min(.data[[v]]), 3),
      .groups = "drop"
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)
}
  
# One-way ANOVA / Kruskal-Wallis to see if sig difference between any start values
for (v in start_vars) {
  kt <- kruskal.test(reformulate("diet", v), data = tank_df)
  cat(sprintf("  %-20s  chi² = %5.2f,  df = %d,  p = %.4f  %s\n",
              v,
              kt$statistic,
              kt$parameter,
              kt$p.value,
              ifelse(kt$p.value < 0.05, "<-- potentially imbalanced", "")))
}

## CHECK ALL STARTING CONDITIONS

# Collapse to one row per tank — start values are constant within tank
tank_df <- df |>
  group_by(tank, diet) |>
  summarise(
    start_ABL     = first(start_ABL),
    start_ABW     = first(start_ABW),
    start_density = first(start_density),
    start_biomass = first(start_biomass),
    .groups = "drop"
  )
 
# Print tank-level table
cat("\n  Starting values per tank:\n")
tank_df |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)
 
# Summary by diet
cat("\n  Starting values summarised by diet:\n")
start_vars <- c("start_ABL", "start_ABW", "start_density", "start_biomass")
for (v in start_vars) {
  cat(sprintf("\n--- %s ---\n", v))
  tank_df |>
    group_by(diet) |>
    summarise(
      n      = n(),
      mean   = round(mean(.data[[v]]),   3),
      sd     = round(sd(.data[[v]]),     3),
      min    = round(min(.data[[v]]),    3),
      max    = round(max(.data[[v]]),    3),
      range  = round(max(.data[[v]]) - min(.data[[v]]), 3),
      .groups = "drop"
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)
}
  
# One-way ANOVA / Kruskal-Wallis to see if sig difference between any start values
for (v in start_vars) {
  kt <- kruskal.test(reformulate("diet", v), data = tank_df)
  cat(sprintf("  %-20s  chi² = %5.2f,  df = %d,  p = %.4f  %s\n",
              v,
              kt$statistic,
              kt$parameter,
              kt$p.value,
              ifelse(kt$p.value < 0.05, "<-- potentially imbalanced", "")))
}

### CORRELATION BETWEEN MODERATORS

## Should include these as moderators in model, check collinearity between them

# Pearson correlation among moderators
mod_mat <- tank_df |>
  dplyr::select(start_ABL, start_ABW, start_density, start_biomass) |>
  as.data.frame()
 
cor_r <- cor(mod_mat, method = "pearson")
cor_p <- cor_pmat(mod_mat)
 
print(round(cor_r, 3))
print(format(round(cor_p, 4)))
 
# Correlation heatmap
p_cor <- ggcorrplot(
  cor_r,
  method    = "square",
  type      = "lower",
  lab       = TRUE,
  lab_size  = 5,
  colors    = c("#C0392B", "white", "#1F618D"),
  title     = "Pearson correlation: tank-level moderators\n(n = 12 tanks)",
  ggtheme   = theme_minimal(base_size = 13)
)
 
print(p_cor)

## start_ABL and start_ABW are correlated, so are start_biomass and density - This makes sense because they are calculated from each other. Proceed with one from each for modelling (start_ABW and density?)

### BAYESION MODELLING
