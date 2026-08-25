setwd("C:/Users/RebeccaPedler/Documents/Ulva_diets_abalone")
library(here)

# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 1: Running Bayesion models

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

### MODEL FITTING CONTROL

## Set refit = TRUE to refit all models from scratch and overwrite saved files
## Set refit = FALSE to load previously fitted models from disk
## Change to TRUE any time you modify model parameters or data

refit <- FALSE

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

## E01 had really high mortality due to TGP spike in late January - remove this tank because of TGP confound

# Remove E01 for all downstream analysis

df_raw <- df_raw |>
  filter(tank != "E01")

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

## Outliers relatively even accross tanks and groups - slightly higher in control tanks

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

# Drop wakame and plot just Ulva vs control
df_no_wakame <- df |>
  filter(diet != "wakame") |>
  mutate(diet = droplevels(diet))

p_wt_log_hist_Ulva <- hist_diet("log_weight", "log(weight) (g)", data = df_no_wakame)
p_wt_log_hist_Ulva

ggsave(here("figures", "p_wt_log_hist_Ulva.png"),plot = p_wt_log_hist_Ulva, dpi = 300, width = 9, height = 8, units = "in")

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

## Animals in the wakame treatment died more, finished at a lower density and had higher feed availability compared to Ulva and commercial
# Any growth improvement cannot be reliably attributed to wakame

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
## end_density and mortality_p are also correlated - again makes sense - less abalone = less density within tanks
## end_density and per_capita_feed weakly correlated 

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

## Tank membership explains approx. 10% of variation between abalone weights. Remaining 90% is individual-level

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

## CREATE FUNCTION FOR PULLING POSTERIOR PROBABILITIES

posterior_summary_table <- function(model, direction = c("greater", "less"), digits = 3) {

  direction <- match.arg(direction)

  draws <- as_draws_df(model) %>%
    select(starts_with("b_"))

  purrr::map_dfr(names(draws), function(param) {
    b <- draws[[param]]

    tibble(
      effect    = param,
      estimate  = median(b),
      est_error = sd(b),
      lower95   = quantile(b, 0.025),
      upper95   = quantile(b, 0.975),
      p_direction = if (direction == "greater") mean(b > 0) else mean(b < 0)
    )
  }) %>%
    rename(!!paste0("p_", direction, "_0") := p_direction) %>%
    mutate(across(where(is.numeric), ~ round(.x, digits)))
}

## SET PRIORS  
 
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

if (refit) {
  fit_weight_default <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z,
    data     = tank_df,
    family   = gaussian(),
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_weight_default, here("models", "fit_weight_default.rds"))
} else {
  fit_weight_default <- readRDS(here("models", "fit_weight_default.rds"))
}

summary(fit_weight_default) 
p_pp_default <- pp_check(fit_weight_default)
ggsave(here("figures", "pp_check_default.png"), plot = p_pp_default, dpi = 300, width = 8, height = 6, units = "in")

# MODEL B2: Gaussian with weakly informative priors

if (refit) {
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
  saveRDS(fit_weight_informative, here("models", "fit_weight_informative.rds"))
} else {
  fit_weight_informative <- readRDS(here("models", "fit_weight_informative.rds"))
}

summary(fit_weight_informative)
p_pp_informative <- pp_check(fit_weight_informative)
ggsave(here("figures", "pp_check_informative.png"), plot = p_pp_informative, dpi = 300, width = 8, height = 6, units = "in")

# MODEL B3: Gaussian with less informative priors (wide normal distribution)

if (refit) {
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
  saveRDS(fit_weight_informative_wide, here("models", "fit_weight_informative_wide.rds"))
} else {
  fit_weight_informative_wide <- readRDS(here("models", "fit_weight_informative_wide.rds"))
}

summary(fit_weight_informative_wide)
p_pp_informative_wide <- pp_check(fit_weight_informative_wide)
ggsave(here("figures", "pp_check_informative_wide.png"), plot = p_pp_informative_wide, dpi = 300, width = 8, height = 6, units = "in")

# MODEL B4: Sensitivity test (add mortality_p) - with weakly informative priors

if (refit) {
  fit_weight_sens <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z + mortality_z,
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
  saveRDS(fit_weight_sens, here("models", "fit_weight_sens.rds"))
} else {
  fit_weight_sens <- readRDS(here("models", "fit_weight_sens.rds"))
}

summary(fit_weight_sens)
p_pp_sens <- pp_check(fit_weight_sens)
ggsave(here("figures", "pp_check_sens.png"), plot = p_pp_sens, dpi = 300, width = 8, height = 6, units = "in")

# Get probabilities 
posterior_summary_table(fit_weight_sens)

# MODEL B5: Sensitivity test (add start_ABW_z) - with weakly informative priors

if (refit) {
  fit_weight_sens2 <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z + start_ABW_z,
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
  saveRDS(fit_weight_sens2, here("models", "fit_weight_sens2.rds"))
} else {
  fit_weight_sens2 <- readRDS(here("models", "fit_weight_sens2.rds"))
}

summary(fit_weight_sens2)
p_pp_sens2 <- pp_check(fit_weight_sens2)
ggsave(here("figures", "pp_check_sens2.png"), plot = p_pp_sens2, dpi = 300, width = 8, height = 6, units = "in")

# Get probabilities 
posterior_summary_table(fit_weight_sens2)

# LEAVE-ONE-OUT COMPARISON

fit_weight_default          <- add_criterion(fit_weight_default,          "loo", moment_match = TRUE)
fit_weight_informative      <- add_criterion(fit_weight_informative,      "loo", moment_match = TRUE)
fit_weight_informative_wide <- add_criterion(fit_weight_informative_wide, "loo", moment_match = TRUE)
fit_weight_sens             <- add_criterion(fit_weight_sens,             "loo", moment_match = TRUE)
fit_weight_sens2            <- add_criterion(fit_weight_sens2,            "loo", moment_match = TRUE)

## Rank all models based on elpd_diff
loo_compare(
  fit_weight_default,
  fit_weight_informative,
  fit_weight_informative_wide,
  fit_weight_sens,
  fit_weight_sens2
)

# Results are robust accross sensitivity analyses and regardless of prior choice

## Continue with model including start_ABW as the final and with more iterations

if (refit) {
  fit_weight_final <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z + start_ABW_z, 
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
  saveRDS(fit_weight_final, here("models", "fit_weight_final.rds"))
} else {
  fit_weight_final <- readRDS(here("models", "fit_weight_final.rds"))
}

summary(fit_weight_final)

# Get probabilities 
posterior_summary_table(fit_weight_final)

### PLOT FINAL MODEL
draws_diet <- fit_weight_final |>
  spread_draws(b_dietulva, b_dietwakame) |>
  mutate(
    ulva   = (exp(b_dietulva)   - 1) * 100,
    wakame = (exp(b_dietwakame) - 1) * 100
  ) |>
  select(.draw, ulva, wakame) |>
  pivot_longer(c(ulva, wakame), names_to = "diet", values_to = "pct_growth") |>
  mutate(diet = factor(diet, levels = c("ulva", "wakame")))

p_effects <- ggplot(draws_diet, aes(x = pct_growth, y = diet, fill = diet)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
  stat_halfeye(
    .width        = c(0.66, 0.95),
    point_interval = median_qi,
    slab_alpha    = 0.75,
    height        = 0.7
  ) +
  scale_fill_manual(values = diet_cols, guide = "none") +
  scale_y_discrete(
    labels = function(x) parse(text = ifelse(x == "ulva", "italic('Ulva')", "'Wakame'"))
  ) +
  labs(
    x = "Increase in final weight relative to control diet (%)",
    y = NULL,
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6)
  )

p_effects
ggsave(here("figures", "p_effects.png"), plot = p_effects, dpi = 300, width = 9, height = 6, units = "in")

## Check diagnostics

# PP check - density overlay
p_pp_density <- pp_check(fit_weight_final)
p_pp_density
ggsave(here("figures", "pp_check_weight_density.png"), plot = p_pp_density, dpi = 300, width = 8, height = 6, units = "in")

# PP check - mean 
p_pp_mean <- pp_check(fit_weight_final, type = "stat", stat = "mean")
p_pp_mean
ggsave(here("figures", "pp_check_weight_mean.png"), plot = p_pp_mean, dpi = 300, width = 6, height = 5, units = "in")

# PP check - SD 
p_pp_sd <- pp_check(fit_weight_final, type = "stat", stat = "sd")
p_pp_sd
ggsave(here("figures", "pp_check_weight_sd.png"), plot = p_pp_sd, dpi = 300, width = 6, height = 5, units = "in")

# Trace plots
p_trace <- mcmc_plot(fit_weight_final, type = "trace")
p_trace
ggsave(here("figures", "trace_weight_final.png"), plot = p_trace, dpi = 300, width = 12, height = 8, units = "in")

# MCMC plots
final_mcmc_trace <- mcmc_trace(fit_weight_final, pars = c("b_Intercept", "b_dietulva", "b_dietwakame", "b_per_capita_feed_z", "sigma"))
final_mcmc_trace
ggsave(here("figures", "final_mcmc_trace.png"), plot = final_mcmc_trace, dpi = 300, width = 12, height = 8, units = "in")

# LOO summary
loo_weight_final <- loo(fit_weight_final)
print(loo_weight_final)

# Check any tanks with pareto_k > 0.7
pareto_k <- loo_weight_final$diagnostics$pareto_k
pareto_k
which(pareto_k > 0.7)

# Do additional sensitivity test for priors

if (refit) {
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
  saveRDS(fit_prior_only, here("models", "fit_prior_only.rds"))
} else {
  fit_prior_only <- readRDS(here("models", "fit_prior_only.rds"))
}

p_pp_prior_only <- pp_check(fit_prior_only) # Confirms that weakly informed prior improve fit
ggsave(here("figures", "pp_check_prior_only.png"), plot = p_pp_prior_only, dpi = 300, width = 8, height = 6, units = "in")

### SUMMARY TABLE

# Control weight
ref_weight <- tank_df |>
  filter(diet == "control") |>
  summarise(ref = mean(exp(mean_log_weight), na.rm = TRUE)) |>
  pull(ref)

b <- as_draws_df(fit_weight_final)$b_dietulva
pct <- exp(b) - 1                          # proportional effect (log weight)
gms <- if (!is.na(ref_weight)) pct * ref_weight else NULL   # absolute grams gain

final_summary_weight <- tibble(
  response    = "log_weight",
  scale       = "log",
  mean        = mean(b),
  median      = median(b),
  sd          = sd(b),
  lower95     = quantile(b, 0.025),
  upper95     = quantile(b, 0.975),
  p_gt_0      = mean(b > 0),
  pct_mean    = mean(pct),
  pct_median  = median(pct),
  pct_lower95 = quantile(pct, 0.025),
  pct_upper95 = quantile(pct, 0.975),
  # absolute grams gain (needs ref_weight)
  g_mean      = if (!is.null(gms)) mean(gms)            else NA_real_,
  g_median    = if (!is.null(gms)) median(gms)          else NA_real_,
  g_lower95   = if (!is.null(gms)) quantile(gms, 0.025) else NA_real_,
  g_upper95   = if (!is.null(gms)) quantile(gms, 0.975) else NA_real_
)

print(final_summary_weight)

## Sensitivity check one - Does the diet effect hold with E04 and E13 removed
# E04 had highest start_ABW
# E13 had the highest start count and lowest per_capita_feed
# Check that Ulva effect is not driven by leverage of covariates

# Create dataset with E04 and E13 removed
sensitivity <- tank_df |>
  filter(!tank %in% c("E04", "E13"))

# Run model

if (refit) {
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
  saveRDS(fit_weight_noinfl, here("models", "fit_weight_noinfl.rds"))
} else {
  fit_weight_noinfl <- readRDS(here("models", "fit_weight_noinfl.rds"))
}

summary(fit_weight_noinfl)
p_pp_noinfl <- pp_check(fit_weight_noinfl)
ggsave(here("figures", "pp_check_noinfl.png"), plot = p_pp_noinfl, dpi = 300, width = 8, height = 6, units = "in")

## Compare the diet effect: full model vs influential-tanks-removed
post_full   <- as_draws_df(fit_weight_final)$b_dietulva
post_noinfl <- as_draws_df(fit_weight_noinfl)$b_dietulva

compare_tbl <- tibble(
  model      = c("Full (11 tanks)", "Excl. E04 + E13 (9 tanks)"),
  median_log = c(median(post_full),  median(post_noinfl)),
  lower95    = c(quantile(post_full, .025),  quantile(post_noinfl, .025)),
  upper95    = c(quantile(post_full, .975),  quantile(post_noinfl, .975)),
  pct_median = c(100 * (exp(median(post_full))   - 1),
                 100 * (exp(median(post_noinfl)) - 1)),
  p_gt_0     = c(mean(post_full > 0), mean(post_noinfl > 0))
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(as.data.frame(compare_tbl), row.names = FALSE)

# Ulva effect unchanged - growth benefit holds up

## Additional sensitivity test 
# Individual-level heirarchical model with tank random effect - does Ulva effect hold up?

# Create dataframe with scaled per_capita_feed
df_sub <- df_sub |>
  mutate(per_capita_feed_z = scale(per_capita_feed)[, 1])

# Run model 

if (refit) {
  fit_hier <- brm(
    log_weight ~ diet + per_capita_feed_z + (1 | tank),
    data    = df_sub,
    family  = gaussian(),
    prior   = priors_logwt,         
    chains  = 4, cores = 4,
    iter    = 4000, warmup = 2000,
    seed    = 42,
    control = list(adapt_delta = 0.95)
  )
  saveRDS(fit_hier, here("models", "fit_hier.rds"))
} else {
  fit_hier <- readRDS(here("models", "fit_hier.rds"))
}

summary(fit_hier)

# Diagnostic checks 
fit_hier_trace <- mcmc_trace(fit_hier, pars = c("b_dietulva", "b_dietwakame", "b_per_capita_feed_z", "sd_tank__Intercept", "sigma"))
ggsave(here("figures", "fit_hier_trace.png"), plot = fit_hier_trace, dpi = 300, width = 6, height = 5, units = "in")
 
fit_hier_ppcheck <- pp_check(fit_hier)
ggsave(here("figures", "fit_hier.png"), plot = fit_hier_ppcheck, dpi = 300, width = 6, height = 5, units = "in")

fit_hier_mean <- pp_check(fit_hier, type = "stat", stat = "mean")
ggsave(here("figures", "fit_hier_mean.png"), plot = fit_hier_mean, dpi = 300, width = 6, height = 5, units = "in")

b_hier <- as_draws_df(fit_hier)$b_dietulva
b_agg  <- as_draws_df(fit_weight_final)$b_dietulva

compare_hier <- tibble(
  model      = c("Aggregated (tank means, n = 12)", "Hierarchical (individuals + 1|tank)"),
  median_log = c(median(b_agg), median(b_hier)),
  lower95    = c(quantile(b_agg, .025), quantile(b_hier, .025)),
  upper95    = c(quantile(b_agg, .975), quantile(b_hier, .975)),
  pct_median = c(100 * (exp(median(b_agg))  - 1), 100 * (exp(median(b_hier)) - 1)),
  p_gt_0     = c(mean(b_agg > 0), mean(b_hier > 0))
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(as.data.frame(compare_hier), row.names = FALSE)

## Check if 10% tank variation from earlier ICC holds
draws    <- as_draws_df(fit_hier)
icc_hier <- draws$sd_tank__Intercept^2 / (draws$sd_tank__Intercept^2 + draws$sigma^2)
c(median = median(icc_hier),
  lower  = quantile(icc_hier, .025),
  upper  = quantile(icc_hier, .975))

# Fit additional model not accounting for per_capita_feed (e.g. simulating the optimisitic growth ceiling for economic analysis if 
# set feed rates had been applied to tanks
## Continue with model including start_ABW as the final and with more iterations

if (refit) {
  fit_weight_final_unadjusted <- brm(
    formula  = mean_log_weight ~ diet + start_ABW_z, 
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
  saveRDS(fit_weight_final_unadjusted, here("models", "fit_weight_final_unadjusted.rds"))
} else {
  fit_weight_final_unadjusted <- readRDS(here("models", "fit_weight_final_unadjusted.rds"))
}

summary(fit_weight_final_unadjusted)

## NEXT: Run economic_modelling.R

### END OF SCRIPT ###
