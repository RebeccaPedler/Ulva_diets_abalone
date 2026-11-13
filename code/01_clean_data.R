# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 1: Inspecting and cleaning data

## Load required packages

### Install packages
install.packages(c("tidyverse","brms","ggplot2","patchwork", "bayesplot","tidybayes","posterior","ggcorrplot", "here", "FSA"))
 
library(tidyverse)
library(brms)
library(ggplot2)
library(patchwork)
library(bayesplot)
library(tidybayes)
library(ggcorrplot)
library(posterior)
library(here)
library(FSA)
library(purrr)
library(dplyr)
library(ggtext)
library(scales)

### MODEL FITTING CONTROL

## Set refit = TRUE to refit all models from scratch and overwrite saved files
## Set refit = FALSE to load previously fitted models from disk
## Change to TRUE any time you modify model parameters or data

refit <- TRUE

### LOAD DATA

df_raw     <- read.csv(here("data", "individual_abalone_data.csv")) # Abalone growth data
df_counts  <- read.csv(here("data","per_capita_feed.csv")) # Running mortality, tank counts, and feed per capita data
str(df_raw)
str(df_counts)

# Make tank and diet factors
df_raw <- df_raw |>
    mutate(
    tank = factor(tank),
    diet = factor(diet, levels = c("control", "ulva", "wakame")),
    log_weight = log(weight_g),
    log_length = log(length_mm)
  )

# Create colour palette for plots
diet_cols <- c(
  "control" = "#8DB4C8",   
  "ulva"    = "#6DAA6E",
  "wakame"  = "#C4845A"
)

### CLEAN AND INSPECT DATA

# Check for missing data
missing_summary <- df_raw |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  mutate(pct_missing = round(n_missing / nrow(df_raw) * 100, 2)) |>
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

# Check sample sizes
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

# Inspect outliers
out_summary <- df |>
  group_by(tank, diet) |>
  mutate(
    z_within = (weight_g - mean(weight_g, na.rm = TRUE)) /
      sd(weight_g, na.rm = TRUE)
  ) |>
  summarise(
    n = sum(!is.na(weight_g)),
    n_outliers = sum(abs(z_within) > 3, na.rm = TRUE),
    percentage_outliers = round((n_outliers / n) * 100, 2),
    .groups = "drop"
  )

print(out_summary)

## Outliers relatively even accross tanks and groups and very small proportion

### CREATE HISTOGRAMS

# Function: histogram + density coloured by diet
hist_diet <- function(var, xlab, binwidth = NULL, data = df) {
  ggplot(data, aes(x = .data[[var]], fill = diet, colour = diet)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins      = 45,
      alpha     = 0.45,
      position  = "identity",
      linewidth = 0.2
    ) +
    geom_density(alpha = 0, linewidth = 0.7) +
    scale_fill_manual(values = diet_cols,   name = "Diet") +
    scale_colour_manual(values = diet_cols, name = "Diet") +
    labs(x = xlab, y = "Density") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      legend.position  = "right"
    )
}

# Function: boxplot by tank, coloured by diet
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
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      legend.position  = "right"
    )
} 

# Print histograms for weight (raw and log)
p_wt_hist <- hist_diet("weight_g", "Final weight (g)") 
p_wt_hist

ggsave(here("figures", "p_wt_hist.png"), plot = p_wt_hist, dpi = 300, width = 9, height = 8, units = "in")
 
p_wt_log_hist <- hist_diet("log_weight", "log(weight) (g)") 
p_wt_log_hist

ggsave(here("figures", "p_wt_log_hist.png"), plot = p_wt_log_hist, dpi = 300, width = 9, height = 8, units = "in")

# Combine plots
p_combined <- p_wt_hist + p_wt_log_hist +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))

p_combined

ggsave(here("figures", "p_wt_combined_hist.png"), plot = p_combined, dpi = 300, width = 14, height = 8, units = "in")

p_wt_box  <- box_tank("weight_g", "Final weight (g)") 
p_wt_box

ggsave(here("figures", "p_wt_box.png"), plot = p_wt_box, dpi = 300, width = 9, height = 8, units = "in")

## Numerical summarry

# Create table
for (v in c("weight_g")) {
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

### COMPARE TANK START AND END CONDITIONS

# Collapse to one row per tank 
tank_df <- df |>
  group_by(tank, diet) |>
  summarise(
    start_density   = first(start_density),
    end_density     = last(end_density),
    start_biomass   = first(start_biomass),
    end_biomass     = first(end_biomass),
    start_count     = first(start_count),
    end_count       = last(end_count),
    start_ABW       = first(start_ABW),
    start_ABL       = first(start_ABL),
    mortality_p     = last(mortality_p),
    per_capita_feed = last(per_capita_feed),
    diet_cost       = first(diet_cost),
    .groups = "drop"
  )
 
# Print tank-level table
tank_df |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Summary by diet
start_vars <- c("start_ABL", "start_ABW", "start_density", "start_biomass", "start_count", 
                "end_density", "end_count", "mortality_p", "per_capita_feed")
for (v in start_vars) {
  cat(sprintf("\n--- %s ---\n", v))
  df |>
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
  cat(sprintf("\n%-20s  chi2 = %5.2f,  df = %d,  p = %.4f  %s\n",
              v,
              as.numeric(kt$statistic),
              as.integer(kt$parameter),
              kt$p.value,
              ifelse(kt$p.value < 0.05, "check", "")))
  # Post-hoc pairwise comparisons if Kruskal-Wallis is significant
  if (kt$p.value < 0.05) {
    dt <- dunnTest(
      reformulate("diet", v),
      data = tank_df,
      method = "bonferroni"
    )
    print(dt$res)
  }
}

### CORRELATION BETWEEN MODERATORS

# Pearson correlation among moderators
mod_mat <- tank_df |>
  dplyr::select(start_ABL, start_ABW, start_density, start_biomass, end_density, mortality_p, per_capita_feed) |>
  as.data.frame()
 
cor_r <- cor(mod_mat, method = "pearson")
cor_p <- cor_pmat(mod_mat)
 
print(round(cor_r, 3))
print(format(round(cor_p, 4)))
 
# Correlation heatmap
p_cor <- ggcorrplot(
  cor_r,
  method   = "square",
  type     = "lower",
  lab      = TRUE,
  lab_size = 5,
  colors   = c("#C0392B", "white", "#1F618D"),
  title    = ""
) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) 
 
print(p_cor)

ggsave(here("figures", "p_cor.png"), plot = p_cor, dpi = 300, width = 12, height = 6, units  = "in")

## start_ABL and start_ABW are correlated, so are start_biomass and density - This makes sense because they are calculated from each other. Proceed with one from each for modelling (start_ABW and density)
## density, mortality and feed per capita also correlated - again makes sense - less abalone = less density within tanks = higher feed availability

### BAYESION MODELLING

## ICC estimates (to determine identify potential tank clustering and see if aggregating to tank level is supported)

# Stratified subsample — 1,000 individuals per tank
set.seed(42)
df_sub <- df |>
  group_by(tank) |>
  slice_sample(n = 1000) |>
  ungroup()

# Run analysis on sub-sample dataset and for abalone weight (g)

if (refit) {
  icc_logweight <- brm(
    formula = log_weight ~ 1 + (1 | tank),
    data    = df_sub,
    family  = gaussian(),
    prior   = c(
      prior(normal(3.7, 0.5),  class = Intercept),  # log scale: exp(3.7) ≈ 40 g
      prior(exponential(1),    class = sd),
      prior(exponential(1),    class = sigma)
    ),
    chains  = 4,
    cores   = 4,
    iter    = 2000,
    warmup  = 1000,
    seed    = 42,
    control = list(adapt_delta = 0.90)
  )
  saveRDS(icc_logweight, here("models", "icc_logweight.rds"))
} else {
  icc_logweight <- readRDS(here("models", "icc_logweight.rds"))
}
 
summary(icc_logweight)
 
var_components_logwt <- VarCorr(icc_logweight)
print(var_components_logwt)
 
draws_logwt   <- as_draws_df(icc_logweight)
sd_tank_logwt  <- draws_logwt$`sd_tank__Intercept`
sd_resid_logwt <- draws_logwt$sigma
 
icc_logwt <- sd_tank_logwt^2 / (sd_tank_logwt^2 + sd_resid_logwt^2)
 
icc_mcmc_trace <- mcmc_trace(icc_logweight, pars = c("b_Intercept", "sd_tank__Intercept", "sigma"))

ggsave(here("figures", "icc_mcmc_trace.png"), plot = icc_mcmc_trace, dpi = 300, width = 12, height = 6, units  = "in")
 
# SUMMARY TABLE
 
icc_summary <- tibble(
  outcome    = c("log_weight"),
  icc_mean   = c(mean(icc_logwt)),
  icc_median = c(median(icc_logwt)),
  icc_lower  = c(quantile(icc_logwt,  0.025)),
  icc_upper  = c(quantile(icc_logwt,  0.975)),
  sd_tank    = c(mean(sd_tank_logwt)),
  sd_resid   = c(mean(sd_resid_logwt))
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))
 
print(icc_summary)

# Tank means differ by ~ 7% in weight and the tank effect is real (CI = 0.05 - 0.12)
# Sigma = 0.35, so abalone within a tank vary by 35% around the tank mean
# ICC = 0.04 - aggregate to tank-level means

## Create tank level dataset for Bayesian models

# Aggregate data to tank level
tank_df <- df |>
  group_by(tank, diet) |>
  summarise(
    mean_log_weight    = mean(log_weight),
    mean_weight_g      = mean(weight_g),
    mean_length_mm     = mean(length_mm),
    start_ABW          = first(start_ABW),
    start_density      = first(start_density),
    start_biomass      = first(start_biomass),
    start_ABL          = first(start_ABL),
    start_count        = first(start_count),
    per_capita_feed    = first(per_capita_feed),
    mortality          = first(mortality_p),
    end_density        = first(end_density),
    diet_cost          = first(diet_cost),
    n                  = n(),
    .groups = "drop"
  )

str(tank_df)

# Save dataset
write_csv(tank_df, here("data", "tank_data.csv"))

## NEXT run "02_bayesian_mortality_models.R"
