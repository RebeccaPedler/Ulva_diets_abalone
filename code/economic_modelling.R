# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 2: Economic analysis

### LOAD PACKAGES
 
library(tidyverse)
library(brms)
library(here)
library(scales)
library(ggtext)
library(ggplot2)
 
### CREATE THEMES
 
# Colour for the single growth model reported here
model_col <- "#0F6E56"
 
# plot theme
theme_ulva <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line.x       = element_line(color = "black", linewidth = 0.6),
      axis.line.y       = element_line(color = "black", linewidth = 0.6),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(-0.2, "cm"),
      legend.position   = "right"
    )
}
 
### LOAD AND CLEAN DATA
 
df_raw <- read.csv(here("data", "individual_abalone_data.csv"))
 
df_raw <- df_raw |>
  mutate(
    tank       = factor(tank),
    diet       = factor(diet, levels = c("control", "ulva", "wakame")),
    log_weight = log(weight_g)
  )
 
## E01 had really high mortality due to TGP spike in late January - remove this tank because
## of TGP confound (matches the exclusion applied in Bayesian_models.R)
df_raw <- df_raw |>
  filter(tank != "E01")
 
# Filter biologically implausible values (segmentation artefacts)
df_econ <- df_raw |>
  filter(
    weight_g >= 10,
    weight_g <= 150
  )
 
cat(sprintf("Rows removed by filter: %d (%d remaining)\n", nrow(df_raw) - nrow(df_econ), nrow(df_econ)))
 
### TANK-LEVEL AGGREGATION
 
tank_df <- df_econ |>
  group_by(tank, diet) |>
  summarise(
    mean_log_weight = mean(log_weight),
    mean_weight_g   = mean(weight_g),
    start_ABW       = first(start_ABW),
    start_count     = first(start_count),
    per_capita_feed = first(per_capita_feed),
    mortality       = first(mortality_p),
    diet_cost       = first(diet_cost),
    n               = n(),
    .groups = "drop"
  ) |>
  mutate(
    per_capita_feed_z = scale(per_capita_feed)[, 1],
    start_ABW_z       = scale(start_ABW)[, 1]
  )
 
cat("\nTank-level summary (control + ulva only):\n")
tank_df |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)
 
### ECONOMIC INPUTS (all prices in AUD)
 
ulva_diet_price <-  3.78   # $/kg finished Ulva diet
ctrl_diet_price <-  2.47   # $/kg finished control diet
ulva_meal_price <-  8.22   # $/kg Ulva meal 
ulva_inclusion  <-  0.20   # Ulva meal inclusion rate (20%)
TARGET_G        <-  90.0   # target harvest weight of abalone 
DAYS            <-  131    # trial duration (days)
FEED_RATE       <-  0.015  # feed rate accross trial period

## Create range for feed's share of total opex (used in sensitivity analysis in Part D)
feed_share_ref  <- 0.20 # Assumes that on average, feed will account for 0.20 of OPEX
feed_share_grid <- c(0.10, 0.15, 0.20, 0.25, 0.30)
 
## Diet premium: extra cost per kg of feed to use the Ulva diet
diet_premium <- ulva_diet_price - ctrl_diet_price   # $1.31/kg
 
## Back out the cost of the ingredient Ulva meal displaces.
## diet_premium = inclusion × (meal_price − displaced_cost)→ displaced_cost = meal_price − diet_premium / inclusion
displaced_cost <- ulva_meal_price - diet_premium / ulva_inclusion
 
cat(sprintf("\n  Diet premium at current prices:    $%.2f/kg feed\n", diet_premium))
cat(sprintf("  Implied displaced ingredient cost: $%.2f/kg meal\n", displaced_cost))
 
### GROWTH MODEL
 
## Primary model: diet + per_capita_feed_z + start_ABW_z
## Adjusting for per-capita feed removes the feed-availability confound
fit_model <- readRDS(here("models", "fit_weight_final.rds"))
summary(fit_model)
 
### EXTRACT POSTERIORS
 
b_growth <- fit_model |>
  as_draws_df() |>
  pull(b_dietulva)
 
growth_pct <- exp(b_growth) - 1
 
cat(sprintf("\n  Growth (primary model):  median %+.2f%%  P(>0) = %.3f\n",
            median(growth_pct) * 100, mean(b_growth > 0)))
 
### PER-ANIMAL DATA FROM TRIAL
## Common starting weight for SGR analysis (control tank mean)
 
W0_ref   <- tank_df |>
  filter(diet == "control") |>
  summarise(w = mean(start_ABW)) |>
  pull(w)
 
## Control mean final weight and observed SGR from the common baseline
ctrl_Wf  <- tank_df |>
  filter(diet == "control") |>
  summarise(w = mean(mean_weight_g)) |>
  pull(w)
 
ctrl_sgr <- (log(ctrl_Wf) - log(W0_ref)) / DAYS   # per day, natural log scale
 
cat(sprintf("\n  Common starting weight (control grand mean):  %.2f g\n", W0_ref))
cat(sprintf("  Control SGR: %.4f %%/day\n", ctrl_sgr * 100))
 
### PART B: DIET PREMIUM AS A FUNCTION OF ULVA MEAL PRICE
## Question: How does diet premium scale with price of raw Ulva meal?
meal_grid_fw <- seq(0, 60, by = 0.5)
 
premium_from_meal <- function(meal_price) {
  ulva_inclusion * (meal_price - displaced_cost)
}
 
premium_per_dollar <- ulva_inclusion
 
cat(sprintf("\n  Each $1/kg rise in Ulva meal price adds $%.2f/kg to the diet premium\n", premium_per_dollar))
cat(sprintf("  At actual meal price ($%.2f/kg): diet premium = $%.2f/kg feed\n", ulva_meal_price, premium_from_meal(ulva_meal_price)))
 
premium_df <- tibble(meal_price = meal_grid_fw) |>
  mutate(
    diet_premium   = premium_from_meal(meal_price),
    ulva_diet_cost = ctrl_diet_price + diet_premium
  )
 
# Readable subset
premium_df |>
  filter(meal_price %in% seq(0, 60, by = 10)) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  as.data.frame() |>
  print(row.names = FALSE)
 
p_premium <- ggplot(premium_df, aes(x = meal_price, y = diet_premium)) +
  geom_vline(xintercept = ulva_meal_price,
             colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0,
             colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 1, colour = "#0F6E56") +
  annotate("text",
           x = ulva_meal_price + 1, y = min(premium_df$diet_premium),
           label = sprintf("Actual meal\nprice $%.2f/kg", ulva_meal_price),
           hjust = 0, vjust = 0, size = 3, colour = "#993C1D") +
  scale_x_continuous(labels = dollar_format()) +
  scale_y_continuous(labels = dollar_format()) +
  labs(
    x        = "*Ulva* meal price ($/kg)",
    y        = "Diet premium over control ($/kg feed)"
  ) +
  theme_ulva() +
    theme(
    axis.title.x = ggtext::element_markdown())
 
print(p_premium)
ggsave(here("figures", "p_premium.png"), plot = p_premium, dpi = 300, width = 8, height = 5, units = "in")
 
### PART C1: TIME-TO-HARVEST BREAK-EVEN USING SGR (OPEX-WEIGHTED)
 
# Ulva final weight per posterior draw, projected from the common baseline W0_ref
ulva_Wf <- function(b) ctrl_Wf * exp(b)
 
# Ulva SGR per draw (estimated from the trial; independent of any harvest target)
ulva_sgr <- function(b) (log(ulva_Wf(b)) - log(W0_ref)) / DAYS
 
# Days to reach a given target weight, for a given SGR
days_to_target <- function(sgr_val, target_g = TARGET_G) {
  log(target_g / W0_ref) / sgr_val
}
 
# Cumulative feed (kg) to reach a given target weight, for a given SGR
feed_to_target_kg <- function(sgr_val, target_g = TARGET_G) {
  FEED_RATE * (target_g - W0_ref) / 1000 / sgr_val
}
 
# Control feed cost and days to a given target
feed_cost_ctrl_at <- function(target_g = TARGET_G) {
  feed_to_target_kg(ctrl_sgr, target_g) * ctrl_diet_price
}
days_ctrl_at <- function(target_g = TARGET_G) {
  days_to_target(ctrl_sgr, target_g)
}
 
# Non-feed opex accrual rate ($/animal/day), benchmarked off the control tank's feed-cost rate and feed_share
non_feed_opex_per_day_for <- function(feed_share = feed_share_ref, target_g_ref = TARGET_G) {
  non_feed_multiplier <- (1 - feed_share) / feed_share
  (feed_cost_ctrl_at(target_g_ref) / days_ctrl_at(target_g_ref)) * non_feed_multiplier
}
 
# Control total cost (feed + non-feed opex) to reach a given target, at a given feed_share
# Equal to feed_cost_ctrl_at(target_g) / feed_share only when target_g == TARGET_G
total_cost_ctrl_at <- function(target_g = TARGET_G, feed_share = feed_share_ref) {
  feed_cost_ctrl_at(target_g) + non_feed_opex_per_day_for(feed_share) * days_ctrl_at(target_g)
}
 
# Ulva diet price implied by a given Ulva meal price (fixed inclusion rate)
ulva_diet_from_meal <- function(meal_price) {
  ctrl_diet_price + ulva_inclusion * (meal_price - displaced_cost)
}
 
# Ulva total cost (feed + non-feed opex) to reach target
total_cost_ulva_at <- function(b, diet_price, target_g = TARGET_G, feed_share = feed_share_ref) {
  feed_to_target_kg(ulva_sgr(b), target_g) * diet_price +
    non_feed_opex_per_day_for(feed_share) * days_to_target(ulva_sgr(b), target_g)
}
 
# Net saving = control total cost − Ulva total cost
total_cost_saving <- function(b, meal_price, target_g = TARGET_G, feed_share = feed_share_ref) {
  total_cost_ctrl_at(target_g, feed_share) -
    total_cost_ulva_at(b, ulva_diet_from_meal(meal_price), target_g, feed_share)
}
 
# Break-even Ulva diet price: diet price at which total cost saving = 0
be_diet_price <- function(b, target_g = TARGET_G, feed_share = feed_share_ref) {
  (total_cost_ctrl_at(target_g, feed_share) -
     non_feed_opex_per_day_for(feed_share) * days_to_target(ulva_sgr(b), target_g)) /
    feed_to_target_kg(ulva_sgr(b), target_g)
}
 
# Break-even Ulva meal price: backed out from the diet price break-even
be_meal_price <- function(b, target_g = TARGET_G, feed_share = feed_share_ref) {
  displaced_cost + (be_diet_price(b, target_g, feed_share) - ctrl_diet_price) / ulva_inclusion
}
 
## Common starting weight and control SGR from the trial
days_saved <- function(b, target_g = TARGET_G) {
  days_ctrl_at(target_g) - days_to_target(ulva_sgr(b), target_g)
}
 
## Summary table (Part C, at the default 90g target and actual $8.22/kg meal price)
feed_cost_ctrl <- feed_cost_ctrl_at(TARGET_G)
total_cost_ctrl <- total_cost_ctrl_at(TARGET_G)
 
summarise_sgr_economics <- function(b, label) {
  saving <- total_cost_saving(b, ulva_meal_price)
  be_m   <- be_meal_price(b)
  be_d   <- be_diet_price(b)
  ds     <- days_saved(b)
  tibble(
    model                  = label,
    sgr_adv_median_pct_day = round(median(b / DAYS * 100), 5),
    days_to_target_ctrl    = round(days_ctrl_at(TARGET_G), 1),
    days_to_target_ulva    = round(median(days_to_target(ulva_sgr(b))), 1),
    days_saved_median      = round(median(ds), 1),
    days_saved_lo95        = round(quantile(ds, 0.025), 1),
    days_saved_hi95        = round(quantile(ds, 0.975), 1),
    p_days_saved_gt0       = round(mean(ds > 0), 3),
    total_cost_ctrl        = round(total_cost_ctrl, 4),
    total_cost_ulva_median = round(median(total_cost_ulva_at(b, ulva_diet_price)), 4),
    net_saving_median      = round(median(saving), 4),
    net_saving_lo95        = round(quantile(saving, 0.025), 4),
    net_saving_hi95        = round(quantile(saving, 0.975), 4),
    p_saving_gt0           = round(mean(saving > 0), 3),
    be_meal_median         = round(median(be_m), 2),
    be_meal_lo95           = round(quantile(be_m, 0.025), 2),
    be_meal_hi95           = round(quantile(be_m, 0.975), 2),
    be_diet_median         = round(median(be_d), 3),
    be_diet_lo95           = round(quantile(be_d, 0.025), 3),
    be_diet_hi95           = round(quantile(be_d, 0.975), 3)
  )
}

# Print summary
sgr_summary <- summarise_sgr_economics(b_growth, "primary_model")
print(as.data.frame(sgr_summary), row.names = FALSE)
 
## PART C2: PROBABILITY-OF-POSITIVE-SAVING CURVE (across Ulva meal price, 90g target)
 
meal_grid_sgr <- seq(0, 20, by = 0.25)

# Create dataframe for plotting
prob_saving_curve <- function(b, label) {
  tibble(meal_price = meal_grid_sgr) |>
    rowwise() |>
    mutate(
      p_saving_pos  = mean(total_cost_saving(b, meal_price) > 0),
      median_saving = median(total_cost_saving(b, meal_price)),
      model         = label
    ) |>
    ungroup()
}
 
saving_curve_df <- prob_saving_curve(b_growth, "primary_model")

# Plot saving curve
p_saving_curve <- ggplot(saving_curve_df,
                         aes(x = meal_price, y = p_saving_pos)) +
  geom_vline(xintercept = ulva_meal_price,
             colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0.95,
             colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_hline(yintercept = 0.50,
             colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 1, colour = model_col) +
  annotate("text",
           x      = ulva_meal_price,
           y      = 0.12,
           label  = sprintf("Actual meal\nprice $%.2f/kg", ulva_meal_price),
           hjust  = -2,
           vjust  = 0,
           size   = 3,
           colour = "#993C1D") +
  scale_x_continuous(labels = dollar_format(), breaks = seq(0, 20, 2)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    x     = "*Ulva* meal price ($/kg)",
    y     = "P(total cost saving > 0)",
    title = ""
  ) +
  theme_ulva() +
  theme(
    axis.title.x = ggtext::element_markdown()
  )
 
print(p_saving_curve)
ggsave(here("figures", "p_sgr_saving_by_meal_price.png"), plot = p_saving_curve, dpi = 300, width = 9, height = 6, units = "in")

# Create tables with highest meal price at which probability of being profitable = 10% to 90% in 10-point steps, then 95% and 99%
confidence_levels <- c(seq(0.10, 0.90, by = 0.10), 0.95, 0.99)

max_price_at_confidence <- function(level, df = saving_curve_df) {
  vals <- df$meal_price[df$p_saving_pos >= level]
  if (length(vals) == 0) NA_real_ else max(vals)
}

thresh <- tibble(
  confidence_level = confidence_levels,
  max_meal_price    = sapply(confidence_levels, max_price_at_confidence)
) |>
  mutate(
    max_meal_price = round(max_meal_price, 2)
  ) |>
  select(confidence_level, max_meal_price)

print(as.data.frame(thresh), row.names = FALSE)
 
### PART C2: POSTERIOR DISTRIBUTION OF DAYS SAVED
 
days_df <- tibble(days_saved = days_saved(b_growth))
 
p_days <- ggplot(days_df, aes(x = days_saved)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins      = 60,
                 alpha     = 0.5,
                 fill      = model_col,
                 colour    = model_col,
                 linewidth = 0.2) +
  geom_density(alpha = 0, colour = model_col, linewidth = 0.8) +
  geom_vline(
    xintercept = median(days_df$days_saved),
    colour = model_col, linewidth = 0.7, linetype = "dashed"
  ) +
  scale_x_continuous(limits = c(-50, 125), breaks = seq(-50, 125, 25)) +
  labs(
    x        = sprintf("Days saved to reach %.0fg harvest weight", TARGET_G),
    y        = "Density"
  ) +
  theme_ulva()
 
print(p_days)
ggsave(here("figures", "p_days_saved.png"), plot = p_days, dpi = 300, width = 9, height = 6, units = "in")
 
### PART D — SENSITIVITY TO HARVEST TARGET WEIGHT 
harvest_grid <- seq(80, 130, by = 10)
 
## Summary for Part D: break-even meal price and profitability at the actual Ulva meal price ($8.22/kg)
summarise_harvest_sensitivity <- function(b, label) {
  tibble(target_g = harvest_grid) |>
    rowwise() |>
    mutate(
      days_saved_median = round(median(days_saved(b, target_g)), 1),
      total_cost_ctrl   = round(total_cost_ctrl_at(target_g), 4),
      be_meal_median    = round(median(be_meal_price(b, target_g)), 2),
      be_meal_lo95      = round(quantile(be_meal_price(b, target_g), 0.025), 2),
      be_meal_hi95      = round(quantile(be_meal_price(b, target_g), 0.975), 2),
      net_saving_median = round(median(total_cost_saving(b, ulva_meal_price, target_g)), 4),
      p_saving_gt0      = round(mean(total_cost_saving(b, ulva_meal_price, target_g) > 0), 3)
    ) |>
    ungroup() |>
    mutate(model = label)
}
 
harvest_sensitivity_df <- summarise_harvest_sensitivity(b_growth, "primary_model")
 
## Plot: break-even meal price as a function of harvest target weight
p_harvest <- ggplot(harvest_sensitivity_df,
                    aes(x = target_g, y = be_meal_median)) +
  geom_hline(yintercept = ulva_meal_price,
             colour = "#993C1D", linewidth = 0.5) +
  geom_ribbon(aes(ymin = be_meal_lo95, ymax = be_meal_hi95),
              alpha = 0.15, colour = NA, fill = model_col) +
  geom_line(linewidth = 1, colour = model_col) +
  scale_x_continuous(breaks = harvest_grid) +
  scale_y_continuous(labels = dollar_format()) +
  labs(
    x        = "Harvest target weight (g)",
    y        = "Break-even *Ulva* meal price ($/kg)"
  ) +
  theme_ulva() +
  theme(
    axis.title.y = ggtext::element_markdown()
  )
 
print(p_harvest)
ggsave(here("figures", "p_harvest_sensitivity.png"), plot = p_harvest, dpi = 300, width = 9, height = 6, units = "in")
 
### PART D — BREAK-EVEN AND 95%-PROFITABLE PRICE ACROSS FEED OPEX SHARE x HARVEST WEIGHT
 
## Question: how sensitive are the break-even price and the 95% profitable price to feed share of total opex, across a range of harvest weights
 
# Create grid
build_price_tables <- function(b, label) {
 
  be_median <- matrix(
    NA_real_, nrow = length(feed_share_grid), ncol = length(harvest_grid),
    dimnames = list(sprintf("feed_%.0f%%", feed_share_grid * 100),
                    sprintf("%.0fg", harvest_grid))
  )
  be_p95 <- be_median
 
  for (i in seq_along(feed_share_grid)) {
    fs <- feed_share_grid[i]
    for (j in seq_along(harvest_grid)) {
      tg <- harvest_grid[j]
      draws <- be_meal_price(b, target_g = tg, feed_share = fs)
      be_median[i, j] <- round(median(draws), 2)
      be_p95[i, j]    <- round(quantile(draws, 0.05), 2)
    }
  }
 
  list(
    breakeven = as.data.frame(be_median) |> tibble::rownames_to_column("feed_opex_share"),
    p95       = as.data.frame(be_p95)    |> tibble::rownames_to_column("feed_opex_share"),
    model     = label
  )
}
 
tables_primary <- build_price_tables(b_growth, "primary_model")
 
## BREAK EVEN PRICE
print(tables_primary$breakeven, row.names = FALSE)
 
## 95% PROFITABLE PRICE
print(tables_primary$p95, row.names = FALSE)
 
### EXPORT ALL DATA AS CSV
write.csv(premium_df,                here("outputs", "ulva_premium_curve.csv"),            row.names = FALSE)
write.csv(sgr_summary,               here("outputs", "ulva_sgr_economic_summary.csv"),     row.names = FALSE)
write.csv(saving_curve_df,           here("outputs", "ulva_sgr_saving_curve.csv"),         row.names = FALSE)
write.csv(thresh,                    here("outputs", "probability_break_evens.csv"),       row.names = FALSE)
write.csv(harvest_sensitivity_df,    here("outputs", "ulva_harvest_sensitivity.csv"),      row.names = FALSE)
write.csv(tables_primary$breakeven,  here("outputs", "ulva_breakeven_by_opexshare.csv"),   row.names = FALSE)
write.csv(tables_primary$p95,        here("outputs", "ulva_p95_by_opexshare.csv"),         row.names = FALSE)
 
### END OF SCRIPT ###
