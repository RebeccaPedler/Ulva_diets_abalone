# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 4: Bayesian SGR models

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

# Trial length (days) used to calculate SGR (10 Nov 2025 to 20 Mar 2026)
trial_days <- 131

# Calculate SGR (% per day) from tank mean log weight and start ABW
# SGR = 100 * (ln(final weight) - ln(start weight)) / days
tank_df <- tank_df |>
  mutate(
    SGR = 100 * (mean_log_weight - log(start_ABW)) / trial_days
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
priors_sgr <- c(
  prior(student_t(3, 0.6, 0.5), class = Intercept),  # Spread between tanks ~0.06, control ~ 0.06 = 0.05 is weakly informative
  prior(normal(0, 0.25),        class = b),  # Expect diet effects that are +/- 0.2 % per day
  prior(exponential(10),        class = sigma)
)

priors_sgr_wide <- c(
  prior(student_t(3, 0.6, 0.5), class = Intercept),
  prior(normal(0, 1),           class = b),
  prior(exponential(1),         class = sigma)
)

# Individual-level model: individual SGR varies more than tank means (SD ~0.2 % per day) - allowing with wider priors
priors_sgr_hier <- c(
  prior(student_t(3, 0.6, 0.5), class = Intercept),
  prior(normal(0, 0.25),        class = b),
  prior(exponential(5),         class = sigma),
  prior(exponential(10),        class = sd)
)

## SGR MODELS

## INITIAL SENSITIVITY CHECK: Does start weight matter (because SGR calculated from each tanks start mean weight)

if (refit) {
  fit_sgr_non_corrected <- brm(
    formula  = SGR ~ diet + per_capita_feed_z, 
    data     = tank_df,
    family   = gaussian(),
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_sgr_non_corrected, here("models", "fit_sgr_non_corrected.rds"))
} else {
  fit_sgr_non_corrected <- readRDS(here("models", "fit_sgr_non_corrected.rds"))
}

summary(fit_sgr_non_corrected)

# MODEL B1: Gaussian with default priors

if (refit) {
  fit_sgr_default <- brm(
    formula  = SGR ~ diet + per_capita_feed_z + start_ABW_z,
    data     = tank_df,
    family   = gaussian(),
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_sgr_default, here("models", "fit_sgr_default.rds"))
} else {
  fit_sgr_default <- readRDS(here("models", "fit_sgr_default.rds"))
}

summary(fit_sgr_default) 
p_pp_default <- pp_check(fit_sgr_default)
ggsave(here("figures", "pp_check_sgr_default.png"), plot = p_pp_default, dpi = 300, width = 8, height = 6, units = "in")

## Start weight and feed availability also matter for SGR (smaller tanks that were fed more have higher SGR)

## Begin modelling with different priors to see how credibility intervals change

## MODEL B2: Gaussian with weakly informative priors

if (refit) {
  fit_sgr_informative <- brm(
    formula  = SGR ~ diet + per_capita_feed_z + start_ABW_z,
    data     = tank_df,
    family   = gaussian(),
    prior    = priors_sgr,
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_sgr_informative, here("models", "fit_sgr_informative.rds"))
} else {
  fit_sgr_informative <- readRDS(here("models", "fit_sgr_informative.rds"))
}

summary(fit_sgr_informative)
p_pp_informative <- pp_check(fit_sgr_informative)
ggsave(here("figures", "pp_check_sgr_informative.png"), plot = p_pp_informative, dpi = 300, width = 8, height = 6, units = "in")

# MODEL B3: Gaussian with less informative priors (wide normal distribution)

if (refit) {
  fit_sgr_informative_wide <- brm(
    formula  = SGR ~ diet + per_capita_feed_z + start_ABW_z,
    data     = tank_df,
    family   = gaussian(),
    prior    = priors_sgr_wide,
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_sgr_informative_wide, here("models", "fit_sgr_informative_wide.rds"))
} else {
  fit_sgr_informative_wide <- readRDS(here("models", "fit_sgr_informative_wide.rds"))
}

summary(fit_sgr_informative_wide)
p_pp_informative_wide <- pp_check(fit_sgr_informative_wide)
ggsave(here("figures", "pp_check_sgr_informative_wide.png"), plot = p_pp_informative_wide, dpi = 300, width = 8, height = 6, units = "in")


# LEAVE-ONE-OUT COMPARISON

fit_sgr_default          <- add_criterion(fit_sgr_default,          "loo", moment_match = TRUE)
fit_sgr_informative      <- add_criterion(fit_sgr_informative,      "loo", moment_match = TRUE)
fit_sgr_informative_wide <- add_criterion(fit_sgr_informative_wide, "loo", moment_match = TRUE)
fit_sgr_non_corrected    <- add_criterion(fit_sgr_non_corrected,    "loo", moment_match = TRUE)

## Rank all models based on elpd_diff
loo_compare(
  fit_sgr_default,
  fit_sgr_informative,
  fit_sgr_informative_wide,
  fit_sgr_non_corrected
)

# fit_SGR_informative has best fit
## Continue with final model (diet + feed availability + start ABW) and with more iterations

if (refit) {
  fit_sgr_final <- brm(
    formula  = SGR ~ diet + per_capita_feed_z + start_ABW_z, 
    data     = tank_df,
    family   = gaussian(),
    prior    = priors_sgr,
    chains   = 4,
    cores    = 4,
    iter     = 4000,
    warmup   = 2000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_sgr_final, here("models", "fit_sgr_final.rds"))
} else {
  fit_sgr_final <- readRDS(here("models", "fit_sgr_final.rds"))
}

summary(fit_sgr_final)

# Get probabilities 
posterior_summary_table(fit_sgr_final)

### PLOT FINAL MODEL
# SGR effects are already on the % per day scale, so no back-transformation is needed
draws_diet <- fit_sgr_final |>
  spread_draws(b_dietulva, b_dietwakame) |>
  mutate(
    ulva   = b_dietulva,
    wakame = b_dietwakame
  ) |>
  select(.draw, ulva, wakame) |>
  pivot_longer(c(ulva, wakame), names_to = "diet", values_to = "sgr_gain") |>
  mutate(diet = factor(diet, levels = c("ulva", "wakame")))

p_effects <- ggplot(draws_diet, aes(x = sgr_gain, y = diet, fill = diet)) +
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
    x = "Increase in SGR relative to control diet (% per day)",
    y = NULL,
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6)
  )

p_effects
ggsave(here("figures", "p_effects_sgr.png"), plot = p_effects, dpi = 300, width = 9, height = 6, units = "in")

## Check diagnostics

# PP check - density overlay
p_pp_density <- pp_check(fit_sgr_final)
p_pp_density
ggsave(here("figures", "pp_check_sgr_density.png"), plot = p_pp_density, dpi = 300, width = 8, height = 6, units = "in")

# PP check - mean 
p_pp_mean <- pp_check(fit_sgr_final, type = "stat", stat = "mean")
p_pp_mean
ggsave(here("figures", "pp_check_sgr_mean.png"), plot = p_pp_mean, dpi = 300, width = 6, height = 5, units = "in")

# PP check - SD 
p_pp_sd <- pp_check(fit_sgr_final, type = "stat", stat = "sd")
p_pp_sd
ggsave(here("figures", "pp_check_sgr_sd.png"), plot = p_pp_sd, dpi = 300, width = 6, height = 5, units = "in")

# Trace plots
p_trace <- mcmc_plot(fit_sgr_final, type = "trace")
p_trace
ggsave(here("figures", "trace_sgr_final.png"), plot = p_trace, dpi = 300, width = 12, height = 8, units = "in")

# MCMC plots
final_mcmc_trace <- mcmc_trace(fit_sgr_final, pars = c("b_Intercept", "b_dietulva", "b_dietwakame", "b_per_capita_feed_z", "b_start_ABW_z", "sigma"))
final_mcmc_trace
ggsave(here("figures", "final_sgr_mcmc_trace.png"), plot = final_mcmc_trace, dpi = 300, width = 12, height = 8, units = "in")

# LOO summary
loo_sgr_final <- loo(fit_sgr_final)
print(loo_sgr_final)

# Check any tanks with pareto_k > 0.7
pareto_k <- loo_sgr_final$diagnostics$pareto_k
pareto_k
which(pareto_k > 0.7)

# Remove each flagged tank individually and then all three
drop_list <- list(
  "all 12 tanks"     = character(0),
  "without C01"      = "C01",
  "without C04"      = "C04",
  "without U04"      = "U04",
  "without all three" = c("C01", "C04", "U04")
)

# Refit the final model for each scenario
fits <- imap(drop_list, function(drop, label) {
  update(
    fit_sgr_final,
    newdata = filter(tank_df, !tank %in% drop),
    chains = 4, iter = 2000, cores = 4, seed = 123
  )
})

# Compare the diet effects
results <- imap_dfr(fits, function(fit, label) {
  draws <- as_draws_df(fit)
  tibble(
    scenario   = label,
    n_tanks    = nrow(fit$data),
    ulva_med   = median(draws$b_dietulva),
    ulva_lower = quantile(draws$b_dietulva, 0.025),
    ulva_upper = quantile(draws$b_dietulva, 0.975),
    p_ulva_pos = mean(draws$b_dietulva > 0),
    wakame_med = median(draws$b_dietwakame),
    wakame_lower = quantile(draws$b_dietwakame, 0.025),
    wakame_upper = quantile(draws$b_dietwakame, 0.975)
  )
})

results # Print

## SGR effect positive in all models but certainty reduces if you remove C01

# Do additional sensitivity test for priors

if (refit) {
  fit_prior_only <- brm(
    formula = SGR ~ diet + per_capita_feed_z + start_ABW_z,
    data    = tank_df,
    family  = gaussian(),
    prior   = priors_sgr,
    sample_prior = "only",
    chains  = 4, 
    cores = 4, 
    iter = 2000, 
    seed = 42
  )
  saveRDS(fit_prior_only, here("models", "fit_sgr_prior_only.rds"))
} else {
  fit_prior_only <- readRDS(here("models", "fit_sgr_prior_only.rds"))
}

p_pp_prior_only <- pp_check(fit_prior_only) # Check that the priors produce plausible SGR values
ggsave(here("figures", "pp_check_sgr_prior_only.png"), plot = p_pp_prior_only, dpi = 300, width = 8, height = 6, units = "in")

### SUMMARY TABLE

# Control SGR (% per day)
ref_sgr <- tank_df |>
  filter(diet == "control") |>
  summarise(ref = mean(SGR, na.rm = TRUE)) |>
  pull(ref)

b <- as_draws_df(fit_sgr_final)$b_dietulva
pct <- if (!is.na(ref_sgr)) 100 * b / ref_sgr else NULL   # SGR gain as % of the control SGR

final_summary_sgr <- tibble(
  response    = "SGR",
  scale       = "% per day",
  mean        = mean(b),
  median      = median(b),
  sd          = sd(b),
  lower95     = quantile(b, 0.025),
  upper95     = quantile(b, 0.975),
  p_gt_0      = mean(b > 0),

  # SGR gain relative to control SGR (needs ref_sgr)
  pct_mean    = if (!is.null(pct)) mean(pct)            else NA_real_,
  pct_median  = if (!is.null(pct)) median(pct)          else NA_real_,
  pct_lower95 = if (!is.null(pct)) quantile(pct, 0.025) else NA_real_,
  pct_upper95 = if (!is.null(pct)) quantile(pct, 0.975) else NA_real_
)

print(final_summary_sgr)

## NEXT: Run "05_economic_modelling.R"

### END OF SCRIPT ###
