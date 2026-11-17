# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 6: Bayesian heavy metal models 
### MODEL FITTING CONTROL

## Set refit = TRUE to refit all models from scratch and overwrite saved files
## Set refit = FALSE to load previously fitted models from disk
## Change to TRUE any time you modify model parameters or data

refit <- TRUE

### LOAD DATA

# Heavy metal concentrations in foot and viscera, one sample of each tissue per tank (12 tanks)
hm <- read.csv(here("data", "heavy_metals.csv"), stringsAsFactors = FALSE)
str(hm)

metals <- c("Arsenic", "Cadmium", "Lead", "Mercury", "Iron")

## Load function for pulling posteriors
source(here("code", "00_helpers.R"))

# Colours match the survival plot in script 04
hm_cols <- c(control = "#1F5FBF", ulva = "#2E9E4F", wakame = "#E8801A")

### PREPARE DATA FOR MODELLING

# Lead values < 0.01, flagged as left-censored e.g, the true value is somewhere below 0.01
hm_long <- hm |>
  mutate(across(all_of(metals), as.character)) |>
  pivot_longer(all_of(metals), names_to = "metal", values_to = "raw") |>
  mutate(
    censored  = if_else(startsWith(raw, "<"), "left", "none"),
    conc      = as.numeric(sub("<", "", raw)),
    log_conc  = log(conc),
    conc_desc = if_else(censored == "left", conc / 2, conc), # half detection limit, descriptive summaries only
    diet      = factor(tolower(diet), levels = c("control", "ulva", "wakame")),
    tissue    = factor(tissue, levels = c("foot", "viscera")),
    tank      = factor(tank),
    metal     = factor(metal, levels = metals)
  )

stopifnot(!anyNA(hm_long$conc))

# Check - how many values are below the detection limit?
hm_long |>
  group_by(metal) |>
  summarise(n = n(), n_censored = sum(censored == "left"), pct_censored = round(100 * mean(censored == "left"), 1))

hm_long |>
  filter(censored == "left") |>
  count(metal, diet, tissue)

### SET PRIORS

int_centres <- hm_long |>
  group_by(metal) |>
  summarise(centre = round(mean(log_conc), 1)) |>
  pull(centre, name = metal)

int_centres

# Weakly informative priors

make_priors <- function(centre, sigma_rate = 2) {
  c(
    set_prior(paste0("student_t(3, ", centre, ", 2)"), class = "Intercept"),
    set_prior("normal(0, 0.5)",                        class = "b"),
    set_prior("normal(0, 2)",                          class = "b", coef = "tissueviscera"),
    set_prior(paste0("exponential(", sigma_rate, ")"), class = "sigma"),
    set_prior(paste0("exponential(", sigma_rate, ")"), class = "sd")
  )
}

make_priors_wide <- function(centre) {
  c(
    set_prior(paste0("student_t(3, ", centre, ", 5)"), class = "Intercept"),
    set_prior("normal(0, 2)",                          class = "b"),
    set_prior("normal(0, 5)",                          class = "b", coef = "tissueviscera"),
    set_prior("exponential(1)",                        class = "sigma"),
    set_prior("exponential(1)",                        class = "sd")
  )
}

sigma_rates <- c(Arsenic = 2, Cadmium = 2, Lead = 1, Mercury = 2, Iron = 2)

priors_hm      <- lapply(metals, function(m) make_priors(int_centres[[m]], sigma_rates[[m]]))
priors_hm_wide <- lapply(metals, function(m) make_priors_wide(int_centres[[m]]))
names(priors_hm)      <- metals
names(priors_hm_wide) <- metals

## MODELS (one model per metal)

f_add <- bf(log_conc | cens(censored) ~ diet + tissue + (1 | tank))
f_int <- bf(log_conc | cens(censored) ~ diet * tissue + (1 | tank))

fits_add  <- list()
fits_int  <- list()

# MODEL M1: Additive (diet + tissue) - diet effect assumed the same in both tissues

for (m in metals) {
  dat <- hm_long |> filter(metal == m)
  file_add <- here("models", paste0("fit_hm_", tolower(m), "_add.rds"))
  
  if (refit) {
    fits_add[[m]] <- brm(
      formula   = log_conc | cens(censored) ~ diet + tissue + (1 | tank),
      data      = dat,
      family    = gaussian(),
      prior     = priors_hm[[m]],
      chains    = 4,
      cores     = 4,
      iter      = 2000,
      warmup    = 1000,
      seed      = 42,
      save_pars = save_pars(all = TRUE),
      control   = list(adapt_delta = 0.95)
    )
    saveRDS(fits_add[[m]], file_add)
  } else {
    fits_add[[m]] <- readRDS(file_add)
  }
}

# MODEL M2: Interaction (diet x tissue) - diet effect allowed to differ between foot and viscera

for (m in metals) {
  dat <- hm_long |> filter(metal == m)
  file_int <- here("models", paste0("fit_hm_", tolower(m), "_int.rds"))
  
  if (refit) {
    fits_int[[m]] <- brm(
      formula   = log_conc | cens(censored) ~ diet * tissue + (1 | tank),
      data      = dat,
      family    = gaussian(),
      prior     = priors_hm[[m]],
      chains    = 4,
      cores     = 4,
      iter      = 2000,
      warmup    = 1000,
      seed      = 42,
      save_pars = save_pars(all = TRUE),
      control   = list(adapt_delta = 0.95)
    )
    saveRDS(fits_int[[m]], file_int)
  } else {
    fits_int[[m]] <- readRDS(file_int)
  }
}

for (m in metals) {
  cat("\n=====", m, "=====\n")
  print(summary(fits_add[[m]]))
  print(summary(fits_int[[m]]))
}

# LEAVE-ONE-OUT COMPARISON 

for (m in metals) {
  fits_add[[m]] <- add_criterion(fits_add[[m]], "loo", moment_match = TRUE)
  fits_int[[m]] <- add_criterion(fits_int[[m]], "loo", moment_match = TRUE)
  
  cat("\n=====", m, "=====\n")
  print(loo_compare(fits_add[[m]], fits_int[[m]]))
}

## Rank models based on elpd_diff 

fits_final <- list()

for (m in metals) {
  file_final <- here("models", paste0("fit_hm_", tolower(m), "_final.rds"))
  
  if (refit) {
    fits_final[[m]] <- update(
      fits_int[[m]],
      chains  = 4,
      cores   = 4,
      iter    = 4000,
      warmup  = 2000,
      seed    = 42,
      control = list(adapt_delta = 0.95) 
    )
    saveRDS(fits_final[[m]], file_final)
  } else {
    fits_final[[m]] <- readRDS(file_final)
  }
}

for (m in metals) {
  cat("\n=====", m, "=====\n")
  print(summary(fits_final[[m]]))
}

## SENSITIVITY CHECK - wider priors (do the diet and tissue effects change?)

fits_wide <- list()

for (m in metals) {
  dat <- hm_long |> filter(metal == m)
  file_wide <- here("models", paste0("fit_hm_", tolower(m), "_wide.rds"))
  
  if (refit) {
    fits_wide[[m]] <- brm(
      formula   = f_int,
      data      = dat,
      family    = gaussian(),
      prior     = priors_hm_wide[[m]],
      chains    = 4,
      cores     = 4,
      iter      = 4000,
      warmup    = 2000,
      seed      = 42,
      save_pars = save_pars(all = TRUE),
      control   = list(adapt_delta = 0.95)
    )
    saveRDS(fits_wide[[m]], file_wide)
  } else {
    fits_wide[[m]] <- readRDS(file_wide)
  }
}

### GET PROBABILITIES

for (m in metals) {
  cat("\n=====", m, "=====\n")
  print(posterior_summary_table(fits_final[[m]]))
}

# Diet effects within each tissue (log scale)
extract_diet_effects <- function(fit, m) {
  d <- as_draws_df(fit)
  tibble(
    metal          = m,
    ulva_foot      = d$b_dietulva,
    wakame_foot    = d$b_dietwakame,
    ulva_viscera   = d$b_dietulva   + d$`b_dietulva:tissueviscera`,
    wakame_viscera = d$b_dietwakame + d$`b_dietwakame:tissueviscera`
  ) |>
    pivot_longer(-metal, names_to = c("diet", "tissue"), names_sep = "_", values_to = "b")
}

# Tissue effect (viscera vs foot) within each diet (log scale)
extract_tissue_effects <- function(fit, m) {
  d <- as_draws_df(fit)
  tibble(
    metal   = m,
    control = d$b_tissueviscera,
    ulva    = d$b_tissueviscera + d$`b_dietulva:tissueviscera`,
    wakame  = d$b_tissueviscera + d$`b_dietwakame:tissueviscera`
  ) |>
    pivot_longer(-metal, names_to = "diet", values_to = "b")
}

diet_effects   <- imap_dfr(fits_final, extract_diet_effects)
tissue_effects <- imap_dfr(fits_final, extract_tissue_effects)

# Diet effect as % change in concentration 
diet_summary <- diet_effects |>
  mutate(pct = (exp(b) - 1) * 100) |>
  group_by(metal, diet, tissue) |>
  summarise(
    median_pct  = median(pct),
    lower95_pct = quantile(pct, 0.025),
    upper95_pct = quantile(pct, 0.975),
    p_gt_0      = mean(b > 0),
    .groups     = "drop"
  ) |>
  mutate(
    metal  = factor(metal, levels = metals),
    diet   = factor(diet, levels = c("ulva", "wakame")),
    tissue = factor(tissue, levels = c("foot", "viscera"))
  ) |>
  arrange(metal, tissue, diet) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(as.data.frame(diet_summary), row.names = FALSE)

# Tissue effect as a fold change 
tissue_summary <- tissue_effects |>
  mutate(fold = exp(b)) |>
  group_by(metal, diet) |>
  summarise(
    median_fold  = median(fold),
    lower95_fold = quantile(fold, 0.025),
    upper95_fold = quantile(fold, 0.975),
    p_gt_0       = mean(b > 0), 
    .groups      = "drop"
  ) |>
  mutate(
    metal = factor(metal, levels = metals),
    diet  = factor(diet, levels = c("control", "ulva", "wakame"))
  ) |>
  arrange(metal, diet) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(as.data.frame(tissue_summary), row.names = FALSE)

### END OF SCRIPT ###
