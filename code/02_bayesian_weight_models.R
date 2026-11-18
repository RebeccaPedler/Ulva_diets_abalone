# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 2: Bayesian weight models

### MODEL FITTING CONTROL

## Set refit = TRUE to refit all models from scratch and overwrite saved files
## Set refit = FALSE to load previously fitted models from disk
## Change to TRUE any time you modify model parameters or data

refit <- TRUE

### LOAD DATA

tank_df     <- read.csv(here("data", "tank_data.csv")) # Abalone growth data
df_counts  <- read.csv(here("data","per_capita_feed.csv")) # Running mortality, tank counts, and feed per capita data
str(tank_df)
str(df_counts)

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

## Load function for pulling posteriors
source(here("code", "00_helpers.R"))

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

## INITIAL SENSITIVITY CHECK: DOES FEED AVAILABILITY MATTER 

if (refit) {
  fit_weight_non_corrected <- brm(
    formula  = mean_log_weight ~ diet, 
    data     = tank_df,
    family   = gaussian(),
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_weight_non_corrected, here("models", "fit_weight_non_corrected.rds"))
} else {
  fit_weight_non_corrected <- readRDS(here("models", "fit_weight_non_corrected.rds"))
}

summary(fit_weight_non_corrected)

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

## Feed availability matters - Ulva effect only shows when accounting for the fact that control tanks were provided with more feed per animal
## Effect of feed behaviour is 6% [4-9%]- Tanks that had less abalone = more feed available = grew more
## Begin modelling with different priors to see how credibility intervals change

## MODEL B2: Gaussian with weakly informative priors

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

# MODEL B4: Sensitivity test (add start_density_z) - with weakly informative priors

if (refit) {
  fit_weight_sens <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z + start_density_z,
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

## Start density has no effect

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
# No effect of start weight

# LEAVE-ONE-OUT COMPARISON

fit_weight_default          <- add_criterion(fit_weight_default,          "loo", moment_match = TRUE)
fit_weight_informative      <- add_criterion(fit_weight_informative,      "loo", moment_match = TRUE)
fit_weight_informative_wide <- add_criterion(fit_weight_informative_wide, "loo", moment_match = TRUE)
fit_weight_non_corrected    <- add_criterion(fit_weight_non_corrected,    "loo", moment_match = TRUE)
fit_weight_sens             <- add_criterion(fit_weight_sens,             "loo", moment_match = TRUE)
fit_weight_sens2            <- add_criterion(fit_weight_sens2,            "loo", moment_match = TRUE)

## Rank all models based on elpd_diff
loo_compare(
  fit_weight_default,
  fit_weight_informative,
  fit_weight_informative_wide,
  fit_weight_non_corrected,
  fit_weight_sens,
  fit_weight_sens2
)

# Results are robust accross sensitivity analyses and regardless of prior choice
## Model that does not include feed availability has the worst fit
## Continue with final model accounting for feed availability and with more iterations

if (refit) {
  fit_weight_final <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z, 
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
fit_weight_final <- readRDS(here("models", "fit_weight_final.rds"))

weight_specs <- c(ulva   = "b_dietulva",
                  wakame = "b_dietwakame",
                  feed   = "b_per_capita_feed_z")

eff_weight <- get_effects(fit_weight_final, weight_specs, transform = pct)

p_weight <- plot_effects(
  eff_weight,
  colours   = c(ulva = diet_cols[["ulva"]], wakame = diet_cols[["wakame"]], feed = "grey60"),
  labels    = c(ulva = "italic('Ulva')", wakame = "'Wakame'", feed = "'Feed availability'"),
  centre    = "median", interval = "qi",      
  xlab      = "Effect on final weight (% change)",
  show_prob = TRUE
  ) +
  labs(tag = "A)") +
  theme(plot.tag = element_text(face = "bold", size = 16))

p_weight
ggsave(here("figures", "p_weight.png"), plot = p_weight, dpi = 300, width = 9, height = 6, units = "in")

## Check diagnostics

# PP check - density overlay
p_pp_density <- pp_check(fit_weight_final)
p_pp_density
ggsave(here("figures", "diagnostics", "pp_check_weight_density.png"), plot = p_pp_density, dpi = 300, width = 8, height = 6, units = "in")

# PP check - mean 
p_pp_mean <- pp_check(fit_weight_final, type = "stat", stat = "mean")
p_pp_mean
ggsave(here("figures", "diagnostics", "pp_check_weight_mean.png"), plot = p_pp_mean, dpi = 300, width = 6, height = 5, units = "in")

# PP check - SD 
p_pp_sd <- pp_check(fit_weight_final, type = "stat", stat = "sd")
p_pp_sd
ggsave(here("figures", "diagnostics", "pp_check_weight_sd.png"), plot = p_pp_sd, dpi = 300, width = 6, height = 5, units = "in")

# MCMC plots
final_mcmc_trace <- mcmc_trace(fit_weight_final, pars = c("b_Intercept", "b_dietulva", "b_dietwakame", "b_per_capita_feed_z", "sigma"))
final_mcmc_trace
ggsave(here("figures", "diagnostics", "final_weight_mcmc_trace.png"), plot = final_mcmc_trace, dpi = 300, width = 12, height = 8, units = "in")

# LOO summary
loo_weight_final <- loo(fit_weight_final)
print(loo_weight_final)

# Check any tanks with pareto_k > 0.7
pareto_k <- loo_weight_final$diagnostics$pareto_k
pareto_k
which(pareto_k > 0.7) # No tanks

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

## SENSITIVITY CHECK - Individual-level heirarchical model with tank random effect - does Ulva effect hold up?

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
fit_hier_ppcheck <- pp_check(fit_hier)
fit_hier_mean <- pp_check(fit_hier, type = "stat", stat = "mean")

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

## NEXT: Run "03_bayesian_SGR_models.R"

### END OF SCRIPT ###
