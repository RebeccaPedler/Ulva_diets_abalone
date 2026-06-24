### Install packages
install.packages(c("tidyverse","brms","ggplot2","patchwork","scales", "bayesplot","tidybayes","posterior","ggcorrplot", "ggpubr"))
 
library(tidyverse)
library(brms)
library(ggplot2)
library(patchwork)
library(bayesplot)
library(tidybayes)
library(ggpubr)
library(ggcorrplot)
library(posterior)
library(scales)

### LOAD DATA

setwd("C:/Users/RebeccaPedler/OneDrive - Yumbah/Documents/R&D/Industry PhD/Trials/Commercial trial/R_datasets")
df_raw     <- read.csv("individual_abalone_data.csv") # Abalone growth data
df_counts  <- read.csv("per_feed_capita.csv") # Running mortality, tank counts, and feed per capita data
str(df_raw)
str(df_counts)

# Make tank and diet factors
df <- df_raw |>
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
  "wakame"  = "#ABA300"
)

### CLEAN AND INSPECT DATA

# Check for missing data
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
 
## Detect outliers and generate numerical summarise

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

### CREATE HISTOGRAMS

# Function: histogram + density coloured by diet
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
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}
 
# length (raw and log)
# Raw
p_len_hist <- hist_diet("length_mm", "Final shell length (mm)") + ggtitle("Shell length distribution by diet")
p_len_hist
 
p_len_box  <- box_tank("length_mm", "Final shell length (mm)") + ggtitle("Shell length by tank (coloured by diet)")
p_len_box

# Log
p_len_log_hist <- hist_diet("log_length", "log(length) (mm)") + ggtitle("log(length) distribution by diet")
p_len_log_hist

# weight (raw and log)
p_wt_hist <- hist_diet("weight_g", "Final weight (g)") + ggtitle("Weight distribution by diet")
p_wt_hist
 
p_wt_log_hist <- hist_diet("log_weight", "log(weight) (g)") + ggtitle("log(Weight) distribution by diet")
p_wt_log_hist
 
p_wt_box  <- box_tank("weight_g", "Final weight (g)") + ggtitle("Weight by tank (coloured by diet)")
p_wt_box

## Numerical summaries

# Define variables
response_vars <- c("length_mm", "weight_g")

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

## CHECK ALL CONDITIONS

# Collapse to one row per tank 
tank_df <- df |>
  group_by(tank, diet) |>
  summarise(
    start_density   = first(start_density),
    end_density     = last(end_density),
    start_biomass   = first(start_biomass),
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
              ifelse(kt$p.value < 0.05, "check", "")))
}

### CORRELATION BETWEEN MODERATORS

## Should include these as moderators in model, check collinearity between them

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
  method    = "square",
  type      = "lower",
  lab       = TRUE,
  lab_size  = 5,
  colors    = c("#C0392B", "white", "#1F618D"),
  title     = "Pearson correlation: tank-level moderators\n(n = 12 tanks)",
  ggtheme   = theme_minimal(base_size = 13)
)
 
print(p_cor)

## start_ABL and start_ABW are correlated, so are start_biomass and density - This makes sense because they are calculated from each other. Proceed with one from each for modelling (start_ABW and density)
## end_density and mortality_p are also correlated - again makes sense - less abalone = less density within tanks
## end_density and per_capita_feed weakly correlated 

### BAYESION MODELLING

## ICC estimates (to determine identify potential tank clustering and see if we can aggregate to tank level)

# Stratified subsample — 1,000 individuals per tank
set.seed(42)
df_sub <- df |>
  group_by(tank) |>
  slice_sample(n = 1000) |>
  ungroup()

cat("Subsampled dataset:", nrow(df_sub), "individuals across", nlevels(df_sub$tank), "tanks\n")

# Run analysis on sub-sample dataset and for abalone length (mm) 

icc <- brm(
  formula = log_weight ~ 1 + (1 | tank),
  data    = df_sub,
  family  = gaussian(),
  prior   = c(
    prior(normal(68, 10),  class = Intercept),
    prior(exponential(1),  class = sd),
    prior(exponential(1),  class = sigma)
  ),
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 42,
  control = list(adapt_delta = 0.90)
)

summary(icc)

var_components <- VarCorr(icc)
print(var_components)

draws_len   <- as_draws_df(icc)
sd_tank_len  <- draws_len$`sd_tank__Intercept`
sd_resid_len <- draws_len$sigma

icc <- sd_tank^2 / (sd_tank^2 + sd_resid^2)

mcmc_trace(icc,
           pars = c("b_Intercept", "sd_tank__Intercept", "sigma"))

# Summary table

icc_summary <- tibble(
  outcome    = c("length_mm"),
  icc_mean   = c(mean(icc_len)),
  icc_median = c(median(icc_len)),
  icc_lower  = c(quantile(icc_len,  0.025)),
  icc_upper  = c(quantile(icc_len,  0.975)),
  sd_tank    = c(mean(sd_tank_len)),
  sd_resid   = c(mean(sd_resid_len))
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(icc_summary)

## Tank membership explains approx. 10% of variation between abalone length

## Bayesian mixed models — aggregated to tank level to conserve proccessing speed
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
    per_capita_feed    = first(per_capita_feed),
    mortality          = first(mortality_p),
    end_density        = first(end_density),
    diet_cost          = first(diet_cost),
    n                  = n(),
    .groups = "drop"
  )

 # Create scaled random effects
tank_df <- tank_df |>
  mutate(
    start_ABW_z         = scale(start_ABW)[, 1],
    start_density_z     = scale(start_density)[, 1],
    start_biomass_z     = scale(start_biomass)[, 1],
    start_ABL_z         = scale(start_ABL)[, 1],
    start_count_z       = scale(start_count)[, 1],
    per_capita_feed_z   = scale(per_capita_feed) [,1],
    mortality_z         = scale(mortality) [,1],
    end_density_z       = scale(end_density)[, 1], 
  )

## CREATE PRIORS FOR ALL MODELS 
 
priors_logwt <- c(
  prior(student_t(3, 3.7, 1), class = Intercept),
  prior(normal(0, 1),         class = b),
  prior(exponential(1),       class = sigma)
)

priors_logwt_wide <- c(
  prior(student_t(3, 3.7, 1), class = Intercept),
  prior(normal(0, 5),         class = b),
  prior(exponential(1),       class = sigma)
)

## WEIGHT MODELS

# MODEL B1: Gaussian with default priors

fit_weight_default <- brm(
  formula  = mean_log_weight ~ diet + per_capita_feed_z,
  data     = tank_df
  family   = gaussian(),
  chains   = 4,
  cores    = 4,
  iter     = 2000,
  warmup   = 1000,
  seed     = 42,
  control  = list(adapt_delta = 0.95)
)

summary(fit_weight_default) 
pp_check(fit_weight_default)

# MODEL B2: Gaussian with weakly informative priors

fit_weight_informative <- brm(
  formula  = mean_log_weight ~ diet + per_capita_feed_z,
  data     = tank_df,
  family   = gaussian(),
  prior    = priors_logwt,
  chains   = 4,
  cores    = 4,
  iter     = 2000,
  warmup   = 1000,
  seed     = 42,
  control  = list(adapt_delta = 0.95)
)

summary(fit_weight_informative)
pp_check(fit_weight_informative)

# MODEL B3: Gaussian with less informative priors (wide normal distribution)

fit_weight_informative_wide <- brm(
  formula  = mean_log_weight ~ diet + per_capita_feed_z,
  data     = tank_df,
  family   = gaussian(),
  prior    = priors_logwt_wide,
  chains   = 4,
  cores    = 4,
  iter     = 2000,
  warmup   = 1000,
  seed     = 42,
  control  = list(adapt_delta = 0.95)
)

summary(fit_weight_informative_wide)
pp_check(fit_weight_informative_wide)

# MODEL B4: Sensitivity test (add mortality_p) - with weakly informative priors
fit_weight_sens <- brm(
  formula  = mean_log_weight ~ diet + per_capita_feed_z + mortality_p_z,
  data     = tank_df,
  family   = gaussian(),
  prior    = priors_logwt,
  chains   = 4,
  cores    = 4,
  iter     = 2000,
  warmup   = 1000,
  seed     = 42,
  control  = list(adapt_delta = 0.95)
)

summary(fit_weight_sens)

# LEAVE-ONE-OUT COMPARISON

fit_weight_default          <- add_criterion(fit_weight_default,          "loo", moment_match = TRUE)
fit_weight_informative      <- add_criterion(fit_weight_informative,      "loo", moment_match = TRUE)
fit_weight_informative_wide <- add_criterion(fit_weight_informative_wide, "loo", moment_match = TRUE)
fit_weight_sens             <- add_criterion(fit_weight_sens,             "loo", moment_match = TRUE)

## Inspect a single model's LOO (elpd_loo, p_loo, Pareto-k diagnostics)
loo(fit_weight_sens)

## Rank all candidates based on elpd_diff
loo_compare(
  fit_weight_default,
  fit_weight_informative,
  fit_weight_informative_wide,
  fit_weight_sens
)

## Continue with model B2 as the final and with more iterations

fit_weight_final <- brm(
  formula  = mean_log_weight ~ diet + per_capita_feed_z, 
  data     = tank_df,
  family   = gaussian(),
  prior    = priors_logwt,
  chains   = 4,
  cores    = 4,
  iter     = 2000,
  warmup   = 1000,
  seed     = 42,
  control  = list(adapt_delta = 0.95)
)

summary(fit_weight_final)

# Check all diagnostics
pp_check(fit_weight_final)
pp_check(fit_weight_final, type = "stat", stat = "mean")
pp_check(fit_weight_final, type = "stat", stat = "sd")
mcmc_plot(fit_weight_final, type = "trace")
loo_weight_final <- loo(fit_weight_final)

# Check any tanks with pareto_k > 0.7
pareto_k <- loo_weight_final$diagnostics$pareto_k
pareto_k
which(pareto_k > 0.7)

# Check tanks E04 and E13
tank_df[c(4, 5), ]

# Refit with moment matching
loo_weight_final_mm <- loo(fit_weight_final, moment_match = TRUE)
loo_weight_final_mm

# Rows 4 and 5 were influential, but they were not so extreme that LOO failed completely. 
# Moment matching adjusted the approximation, and the overall estimates barely changed. That is reassuring.

# Do additional sensitivity test for priors
fit_prior_only <- brm(
  formula = mean_log_weight ~ diet + per_capita_feed_z,
  data    = tank_df,
  family  = gaussian(),
  prior   = priors_logwt,
  sample_prior = "only",
  chains  = 4, 
  cores = 4, 
  iter = 2000, 
  seed = 42
)
pp_check(fit_prior_only) # Confirms that weakly informed prior improve fit

## Summary table for both models
## Summary table for both models — log scale, proportional (%), and grams
summarise_final <- function(fit, response, scale = c("identity", "log"),
                            ref_weight = NA_real_) {
  scale <- match.arg(scale)
  b <- as_draws_df(fit)$b_dietulva

  # proportional effect (log weight)
  pct <- if (scale == "log") exp(b) - 1 else NULL

  # absolute grams gain = proportional effect x reference control weight (g)
  gms <- if (scale == "log" && !is.na(ref_weight)) pct * ref_weight else NULL
  tibble(
    response    = response,
    scale       = scale,
    # --- model (link) scale ---
    mean        = mean(b),
    median      = median(b),
    sd          = sd(b),
    lower95     = quantile(b, 0.025),
    upper95     = quantile(b, 0.975),
    p_gt_0      = mean(b > 0),
    # proportional (%) effect: log model only 
    pct_mean    = if (scale == "log") mean(pct)            else NA_real_,
    pct_median  = if (scale == "log") median(pct)          else NA_real_,
    pct_lower95 = if (scale == "log") quantile(pct, 0.025) else NA_real_,
    pct_upper95 = if (scale == "log") quantile(pct, 0.975) else NA_real_,
    # absolute grams gain: log model only, needs ref_weight 
    g_mean      = if (!is.null(gms)) mean(gms)            else NA_real_,
    g_median    = if (!is.null(gms)) median(gms)          else NA_real_,
    g_lower95   = if (!is.null(gms)) quantile(gms, 0.025) else NA_real_,
    g_upper95   = if (!is.null(gms)) quantile(gms, 0.975) else NA_real_
  )
}

## Reference control weight (g) — the baseline the % gain is applied to.
ref_control_weight <- df |>
  filter(diet == "control") |>
  summarise(w = mean(weight_g)) |>
  pull(w)

final_summary <- bind_rows(
  summarise_final(fit_length_final, "length_mm",  scale = "identity"),
  summarise_final(fit_weight_final, "log_weight", scale = "log",
                  ref_weight = ref_control_weight)
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(as.data.frame(final_summary), row.names = FALSE)

## Sensitivity check one - Does the diet effect hold with E04 and E13 removed

# Create dataset with E04 and E13 removed
sensitivity <- tank_df |>
  filter(!tank %in% c("E04", "E13"))

# Run model
fit_weight_noinfl <- brm(
  formula  = mean_log_weight ~ diet + per_capita_feed_z,
  data     = sensitivity,
  family   = gaussian(),
  prior    = priors_logwt,
  chains   = 4,
  cores    = 4,
  iter     = 4000,
  warmup   = 2000,
  seed     = 42,
  control  = list(adapt_delta = 0.99, max_treedepth = 12)
)

summary(fit_weight_noinfl)
pp_check(fit_weight_noinfl)

## Compare the diet effect: full model vs influential-tanks-removed
post_full   <- as_draws_df(fit_weight_final)$b_dietulva
post_noinfl <- as_draws_df(fit_weight_noinfl)$b_dietulva

compare_tbl <- tibble(
  model      = c("Full (8 tanks)", "Excl. E04 + E13 (6 tanks)"),
  median_log = c(median(post_full),  median(post_noinfl)),
  lower95    = c(quantile(post_full, .025),  quantile(post_noinfl, .025)),
  upper95    = c(quantile(post_full, .975),  quantile(post_noinfl, .975)),
  pct_median = c(100 * (exp(median(post_full))   - 1),
                 100 * (exp(median(post_noinfl)) - 1)),
  p_gt_0     = c(mean(post_full > 0), mean(post_noinfl > 0))
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(as.data.frame(compare_tbl), row.names = FALSE)

### Next - economic break-even

# Fit model for weight without correcting for feed availability
fit_weight_informative_uncorrected <- brm(
  formula  = mean_log_weight ~ diet
  data     = tank_df,
  family   = gaussian(),
  prior    = priors_logwt,
  chains   = 4,
  cores    = 4,
  iter     = 4000,
  warmup   = 2000,
  seed     = 42,
  control  = list(adapt_delta = 0.95)
)

summary(fit_weight_informative_uncorrected)
