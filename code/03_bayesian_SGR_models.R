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
ggsave(here("figures", "pp_check_default.png"), plot = p_pp_default, dpi = 300, width = 8, height = 6, units = "in")

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
ggsave(here("figures", "pp_check_sens.png"), plot = p_pp_sens, dpi = 300, width = 8, height = 6, units = "in")

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
ggsave(here("figures", "pp_check_sens2.png"), plot = p_pp_sens2, dpi = 300, width = 8, height = 6, units = "in")

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

## Sensitivity check one - Does the diet effect hold with E02 removed (Pareto K > 0.7)

# Create dataset with E02 removed
sensitivity <- tank_df |>
  filter(!tank %in% c("E02"))

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
  model      = c("Full (11 tanks)", "Excl. E02 (10 tanks)"),
  median_log = c(median(post_full),  median(post_noinfl)),
  lower95    = c(quantile(post_full, .025),  quantile(post_noinfl, .025)),
  upper95    = c(quantile(post_full, .975),  quantile(post_noinfl, .975)),
  pct_median = c(100 * (exp(median(post_full))   - 1),
                 100 * (exp(median(post_noinfl)) - 1)),
  p_gt_0     = c(mean(post_full > 0), mean(post_noinfl > 0))
) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(as.data.frame(compare_tbl), row.names = FALSE)

# Ulva effect increases slightly - growth benefit holds up

## Sensitivity check two - Individual-level heirarchical model with tank random effect - does Ulva effect hold up?

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


## NEXT: Run "03_bayesian_SGR_models.R"

### END OF SCRIPT ###
