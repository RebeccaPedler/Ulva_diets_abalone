# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 4: Did mortality affect growth? (log final weight and SGR)

### MODEL FITTING CONTROL

## Set refit = TRUE to refit the mortality models from scratch and overwrite saved files
## Set refit = FALSE to load previously fitted mortality models from disk
## The base models are always loaded from disk

refit <- TRUE

### LOAD DATA

tank_df     <- read.csv(here("data", "tank_data.csv")) 
df_counts <- read.csv(here("data", "per_capita_feed.csv"))          
df_ind    <- read.csv(here("data", "individual_abalone_data.csv"))  

df_counts <- df_counts |>
  mutate(date = dmy(date))

str(df_counts)

# Trial length (days) used to calculate SGR (10 Nov 2025 to 20 Mar 2026)
trial_days <- 131

# Calculate SGR 
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

# Mortality by diet
tank_df |>
  group_by(diet) |>
  summarise(mean_mortality = mean(mortality), sd_mortality = sd(mortality), .groups = "drop")

# Correlation between mortality and the other scaled covariates
round(cor(tank_df[, c("mortality_z", "per_capita_feed_z", "start_density_z", "end_density_z", "start_ABW_z")]), 2)

## Load function for pulling posteriors
source(here("code", "00_helpers.R"))

## LOG FINAL WEIGHT MODELS

# BASE MODEL: final weight model (diet + feed availability)
fit_weight_final <- readRDS(here("models", "fit_weight_final.rds"))

summary(fit_weight_final)

# MODEL M1: Diet + mortality (total effect of mortality, feed availability not included)

if (refit) {
  fit_mort_weight_total <- brm(
    formula  = mean_log_weight ~ diet + mortality_z,
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
  saveRDS(fit_mort_weight_total, here("models", "fit_mort_weight_total.rds"))
} else {
  fit_mort_weight_total <- readRDS(here("models", "fit_mort_weight_total.rds"))
}

summary(fit_mort_weight_total)

# MODEL M2: Diet + feed availability + mortality (direct effect of mortality)

if (refit) {
  fit_mort_weight_direct <- brm(
    formula  = mean_log_weight ~ diet + per_capita_feed_z + mortality_z,
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
  saveRDS(fit_mort_weight_direct, here("models", "fit_mort_weight_direct.rds"))
} else {
  fit_mort_weight_direct <- readRDS(here("models", "fit_mort_weight_direct.rds"))
}

summary(fit_mort_weight_direct)

## SGR MODELS

# BASE MODEL: final SGR model (diet + feed availability + start ABW)
fit_sgr_final <- readRDS(here("models", "fit_sgr_final.rds"))

summary(fit_sgr_final)

# MODEL S1: Diet + start ABW + mortality (total effect of mortality, feed availability not included)

if (refit) {
  fit_mort_sgr_total <- brm(
    formula  = SGR ~ diet + mortality_z + start_ABW_z,
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
  saveRDS(fit_mort_sgr_total, here("models", "fit_mort_sgr_total.rds"))
} else {
  fit_mort_sgr_total <- readRDS(here("models", "fit_mort_sgr_total.rds"))
}

summary(fit_mort_sgr_total)

# MODEL S2: Diet + feed availability + start ABW + mortality (direct effect of mortality)

if (refit) {
  fit_mort_sgr_direct <- brm(
    formula  = SGR ~ diet + per_capita_feed_z + start_ABW_z + mortality_z,
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
  saveRDS(fit_mort_sgr_direct, here("models", "fit_mort_sgr_direct.rds"))
} else {
  fit_mort_sgr_direct <- readRDS(here("models", "fit_mort_sgr_direct.rds"))
}

summary(fit_mort_sgr_direct)

### GET PROBABILITIES

# Weight: compare diet effects with and without mortality, and look at the mortality_z row
posterior_summary_table(fit_weight_final)
posterior_summary_table(fit_mort_weight_total)
posterior_summary_table(fit_mort_weight_direct)

# SGR
posterior_summary_table(fit_sgr_final)
posterior_summary_table(fit_mort_sgr_total)
posterior_summary_table(fit_mort_sgr_direct)

### LEAVE-ONE-OUT COMPARISON (does adding mortality improve out-of-sample prediction?)

fit_weight_final   <- add_criterion(fit_weight_final,   "loo", moment_match = TRUE)
fit_mort_weight_total  <- add_criterion(fit_mort_weight_total,  "loo", moment_match = TRUE)
fit_mort_weight_direct <- add_criterion(fit_mort_weight_direct, "loo", moment_match = TRUE)

fit_sgr_final   <- add_criterion(fit_sgr_final,   "loo", moment_match = TRUE)
fit_mort_sgr_total  <- add_criterion(fit_mort_sgr_total,  "loo", moment_match = TRUE)
fit_mort_sgr_direct <- add_criterion(fit_mort_sgr_direct, "loo", moment_match = TRUE)

## Rank all models based on elpd_diff
loo_compare(fit_weight_final, fit_mort_weight_total, fit_mort_weight_direct)
loo_compare(fit_sgr_final, fit_mort_sgr_total, fit_mort_sgr_direct)

## No evidence that mortality influenced growth - best fitting models are still the final models from script 02 and 03

## MORTALITY MODELS - HOW DO DIET, TEMPERATURE AND TANK FACTORS INFLUENCE SUMMER MORTALITY

# Check overlap between diet and start covariates
tank_df |>
  group_by(diet) |>
  summarise(mean_start_density = mean(start_density), mean_start_ABW = mean(start_ABW), .groups = "drop")

### WEEKLY INTERVALS

# Create weekly intervals
start_date <- as.Date("2025-11-10")
end_date   <- as.Date("2026-03-20")
bounds     <- sort(unique(c(seq(start_date, end_date, by = "7 days"), end_date)))

counts_long <- df_counts |>
  select(date, ends_with("_count")) |>
  pivot_longer(-date, names_to = "tank", values_to = "n_alive", names_pattern = "(.*)_count")

# Mean water temperature for each weekly interval (scaled)
env_interval <- df_counts |>
  select(date, min_water_temperature, ave_water_temperature, max_water_temperature) |>
  mutate(interval_start = bounds[pmin(findInterval(as.numeric(date), as.numeric(bounds)), length(bounds) - 1)]) |>
  group_by(interval_start) |>
  summarise(
    temp_cover = mean(!is.na(ave_water_temperature)),
    temp_mean  = mean(ave_water_temperature, na.rm = TRUE),
    temp_min   = mean(min_water_temperature, na.rm = TRUE),
    temp_max   = mean(max_water_temperature, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(temp_mean = if_else(temp_cover >= 0.5, temp_mean, NA_real_))

env_interval <- env_interval |>
  mutate(temp_z = scale(temp_mean)[, 1]) # Scale average water temperature

# Check
summary(env_interval)
str(env_interval)

## Summary water temperature

# Weekly
cat(sprintf("\nWater temperature: weekly mean %.1f C (range %.1f to %.1f); 1 SD = %.2f C\n",
            mean(env_interval$temp_mean), min(env_interval$temp_mean), max(env_interval$temp_mean), sd(env_interval$temp_mean)))

# Whole trial (daily readings)
cat(sprintf("Water temperature: whole trial mean %.1f C (SE %.2f; range %.1f to %.1f; n = %d daily readings, %s to %s)\n",
            mean(df_counts$ave_water_temperature, na.rm = TRUE),
            se(df_counts$ave_water_temperature, na.rm = TRUE),
            min(df_counts$min_water_temperature, na.rm = TRUE),
            max(df_counts$max_water_temperature, na.rm = TRUE),
            sum(!is.na(df_counts$ave_water_temperature)),
            format(min(df_counts$date)), format(max(df_counts$date))))

### PREPARE DATA FOR MODELLING

# Tank-level covariates (scaled across the 12 tanks in the earlier script)
tank_cov <- tank_df |>
  select(tank, diet, start_density_z, start_ABW_z) |>
  mutate(diet = factor(diet, levels = c("control", "ulva", "wakame")))

# Check
stopifnot(nrow(tank_cov) == 12)

# Deaths in each tank and week (at risk = animals alive at the start of the week)
mort <- counts_long |>
  filter(date %in% bounds) |>
  arrange(tank, date) |>
  group_by(tank) |>
  mutate(
    interval_start = lag(date),
    at_risk        = round(lag(n_alive)),
    deaths         = at_risk - round(n_alive),
    days           = as.numeric(date - interval_start)
  ) |>
  ungroup() |>
  filter(!is.na(deaths)) |>
  mutate(
    week  = factor(match(interval_start, bounds)),
    month = factor(month.abb[lubridate::month(interval_start + days / 2)],
                   levels = c("Nov", "Dec", "Jan", "Feb", "Mar"))
  ) |>
  left_join(tank_cov, by = "tank") |>
  left_join(env_interval |> select(interval_start, temp_mean, temp_min, temp_max, temp_z), by = "interval_start") |>
  mutate(
    tank       = factor(tank),
    obs        = factor(row_number()),                       # tank x week, for the observation-level random effect
    hazard_pct = 100 * deaths / at_risk / days               # % of animals alive dying per day
  )

# Check
stopifnot(nrow(mort) == 12 * (length(bounds) - 1))
stopifnot(all(mort$deaths >= 0), all(mort$deaths <= mort$at_risk))

# Deaths by month and diet
mort |>
  group_by(month, diet) |>
  summarise(deaths = sum(deaths), .groups = "drop") |>
  pivot_wider(names_from = diet, values_from = deaths)

# Temperature by month, and how much of the week-to-week temperature variation is explained by month
week_df <- mort |>
  distinct(week, month, temp_mean)

week_df |>
  group_by(month) |>
  summarise(n_weeks = n(), mean_temp = mean(temp_mean), min_temp = min(temp_mean), max_temp = max(temp_mean), .groups = "drop")

cat(sprintf("\nR-squared of weekly mean temperature explained by month: %.2f\n",
            summary(lm(temp_mean ~ month, data = week_df))$r.squared))

# Water temp and month are correlated and mortality tracks water temperature

## SET PRIORS

priors_mort <- c(
  prior(student_t(3, -7, 2.5), class = Intercept),
  prior(normal(0, 1),          class = b),
  prior(exponential(2),        class = sd)
)

## MODELS (all fitted to all weeks)


## STEP A: TEMPERATURE VS MONTH (decide which one goes into the other models)

# MODEL H0_1: Water temperature only
if (refit) {
  fit_mort_H0_1 <- brm(
    formula = deaths | trials(at_risk) ~ offset(log(days)) + temp_z + (1 | tank) + (1 | week) + (1 | obs),
    data = mort, family = binomial(link = "cloglog"), prior = priors_mort,
    chains = 4, cores = 4, iter = 2000, warmup = 1000, seed = 42, control = list(adapt_delta = 0.95))
  saveRDS(fit_mort_H0_1, here("models", "fit_mort_H0_1.rds"))
} else {
  fit_mort_H0_1 <- readRDS(here("models", "fit_mort_H0_1.rds"))
}
summary(fit_mort_H0_1)

# MODEL H0_2: Month only
if (refit) {
  fit_mort_H0_2 <- brm(
    formula = deaths | trials(at_risk) ~ offset(log(days)) + month + (1 | tank) + (1 | week) + (1 | obs),
    data = mort, family = binomial(link = "cloglog"), prior = priors_mort,
    chains = 4, cores = 4, iter = 2000, warmup = 1000, seed = 42, control = list(adapt_delta = 0.95))
  saveRDS(fit_mort_H0_2, here("models", "fit_mort_H0_2.rds"))
} else {
  fit_mort_H0_2 <- readRDS(here("models", "fit_mort_H0_2.rds"))
}
summary(fit_mort_H0_2)

# Compare if water temperature or month describes the weekly pattern of mortality better?
fit_mort_H0_1 <- add_criterion(fit_mort_H0_1, "loo")
fit_mort_H0_2 <- add_criterion(fit_mort_H0_2, "loo")
loo_compare(fit_mort_H0_1, fit_mort_H0_2)

## STEP B: CONTINUE WITH WATER TEMPERATURE 

# MODEL H1: WATER TEMP + diet

if (refit) {
  fit_mort_diet <- brm(
    formula  = deaths | trials(at_risk) ~ offset(log(days)) + temp_z + diet + (1 | tank) + (1 | week) + (1 | obs),
    data     = mort,
    family   = binomial(link = "cloglog"),
    prior    = priors_mort,
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_mort_diet, here("models", "fit_mort_diet.rds"))
} else {
  fit_mort_diet <- readRDS(here("models", "fit_mort_diet.rds"))
}

summary(fit_mort_diet)

# MODEL H2: WATER TEMP + start density

if (refit) {
  fit_mort_density <- brm(
    formula  = deaths | trials(at_risk) ~ offset(log(days)) + temp_z + start_density_z + (1 | tank) + (1 | week) + (1 | obs),
    data     = mort,
    family   = binomial(link = "cloglog"),
    prior    = priors_mort,
    chains   = 4,
    cores    = 4,
    iter     = 2000,
    warmup   = 1000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_mort_density, here("models", "fit_mort_density.rds"))
} else {
  fit_mort_density <- readRDS(here("models", "fit_mort_density.rds"))
}
summary(fit_mort_density)

# MODEL H3: WATER TEMP + start ABW

if (refit) {
  fit_mort_ABW <- brm(
    formula  = deaths | trials(at_risk) ~ offset(log(days)) + month + start_ABW_z + (1 | tank) + (1 | week) + (1 | obs),
    data     = mort,
    family   = binomial(link = "cloglog"),
    prior    = priors_mort,
    chains   = 4,
    cores    = 4,
    iter     = 4000,
    warmup   = 2000,
    seed     = 42,
    control  = list(adapt_delta = 0.99)
  )
  saveRDS(fit_mort_ABW, here("models", "fit_mort_ABW.rds"))
} else {
  fit_mort_ABW <- readRDS(here("models", "fit_mort_ABW.rds"))
}

summary(fit_mort_ABW)

# MODEL H1: WATER TEMP (final with more iterations)

if (refit) {
  fit_mort_final <- brm(
    formula  = deaths | trials(at_risk) ~ offset(log(days)) + temp_z + (1 | tank) + (1 | week) + (1 | obs),
    data     = mort,
    family   = binomial(link = "cloglog"),
    prior    = priors_mort,
    chains   = 4,
    cores    = 4,
    iter     = 4000,
    warmup   = 2000,
    seed     = 42,
    control  = list(adapt_delta = 0.95)
  )
  saveRDS(fit_mort_final, here("models", "fit_mort_final.rds"))
} else {
  fit_mort_final <- readRDS(here("models", "fit_mort_final.rds"))
}

summary(fit_mort_final)

### PLOT MORTALITY

# Survival (% of starting count) for each tank at the end of each week
surv_tank <- mort |>
  arrange(tank, date) |>
  group_by(tank) |>
  mutate(survival = 100 * n_alive / first(at_risk)) |>
  ungroup() |>
  select(tank, diet, date, survival)

surv_start <- mort |>
  group_by(tank, diet) |>
  summarise(date = min(interval_start), survival = 100, .groups = "drop")

surv_diet <- bind_rows(surv_start, surv_tank) |>
  group_by(diet, date) |>
  summarise(value = mean(survival), .groups = "drop") |>
  mutate(series = recode(as.character(diet), control = "Control", ulva = "Ulva", wakame = "Wakame"))

# Daily water temperature straight from the daily readings
temp_daily <- df_counts |>
  transmute(date = as.Date(date),
            temp_mean = ave_water_temperature,
            temp_min  = min_water_temperature,
            temp_max  = max_water_temperature) |>
  filter(!is.na(temp_mean))

# Map temperature onto survival axis
surv_lim <- c(60, 100)
temp_lim <- c(16, 24)    # covers the full trial range (16.7 to 23.9 C)
temp_to_surv <- function(t) surv_lim[1] + (t - temp_lim[1]) / diff(temp_lim) * diff(surv_lim)
surv_to_temp <- function(s) temp_lim[1] + (s - surv_lim[1]) / diff(surv_lim) * diff(temp_lim)

temp_df <- temp_daily |>
  transmute(date,
            value  = temp_to_surv(temp_mean),
            ymin   = temp_to_surv(temp_min),
            ymax   = temp_to_surv(temp_max),
            series = "Water temperature")

plot_df <- bind_rows(surv_diet |> select(date, value, series),
                     temp_df   |> select(date, value, series)) |>
  mutate(series = factor(series, levels = c("Control", "Ulva", "Wakame", "Water temperature")))

p_survival <- ggplot(plot_df, aes(date, value, colour = series, linetype = series)) +
  # daily min-max band drawn first so it sits behind the lines
  geom_ribbon(data = temp_df, aes(x = date, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(Control = "#1F5FBF", Ulva = "#2E9E4F", Wakame = "#E8801A",
                                 `Water temperature` = "black")) +
  scale_linetype_manual(values = c(Control = "solid", Ulva = "dotted", Wakame = "dashed",
                                   `Water temperature` = "solid")) +
  scale_x_date(breaks = seq(as.Date("2025-11-01"), as.Date("2026-03-01"), by = "month"),
               labels = date_format("%b"),
               limits = c(start_date, end_date), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = c(0, 0),
                     sec.axis = sec_axis(~ surv_to_temp(.), name = "Water temperature (°C)")) +
  coord_cartesian(ylim = surv_lim) +
  labs(x = "Month", y = "Survival (%)", colour = NULL, linetype = NULL) +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")

p_survival
ggsave(here("figures", "p_survival_temperature.png"), p_survival, dpi = 300, width = 8, height = 5)

## NEXT: Run "05_economic_modelling.R"

### END OF SCRIPT ###
