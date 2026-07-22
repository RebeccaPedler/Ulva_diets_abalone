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

# Colour palette
model_cols <- c(
  "unadjusted" = "#185FA5",
  "adjusted"   = "#0F6E56"
)

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

# Filter biologically implausible values (segmentation artefacts)
df <- df_raw |>
  filter(
    weight_g >= 10,
    weight_g <= 150
  )

cat(sprintf("Rows removed by filter: %d (%d remaining)\n",
            nrow(df_raw) - nrow(df), nrow(df)))

### TANK-LEVEL AGGREGATION

tank_df <- df_econ |>
  group_by(tank, diet) |>
  summarise(
    mean_log_weight = mean(log_weight),
    mean_weight_g   = mean(weight_g),
    start_ABW       = first(start_ABW),
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

### ECONOMIC INPUTS
## All prices AUD. Feed quantities converted to kg where relevant.

farmgate_price  <- 25.00   # $/kg whole abalone at farmgate
ulva_diet_price <-  3.78   # $/kg finished Ulva diet
ctrl_diet_price <-  2.47   # $/kg finished control diet
ulva_meal_price <-  8.22   # $/kg Ulva meal (the key negotiable lever)
ulva_inclusion  <-  0.20   # Ulva meal inclusion rate (20%)
TARGET_G        <- 80.0    # target harvest weight for SGR analysis (g, whole animal)
DAYS            <- 131     # trial duration (days)
FEED_RATE       <- 0.01    # assumed feed rate for SGR analysis (1% BW/day;
                           # cancels algebraically in net cost comparisons)

## Diet premium: extra cost per kg of feed to use the Ulva diet
diet_premium <- ulva_diet_price - ctrl_diet_price   # $1.31/kg

## Back out the cost of the ingredient Ulva meal displaces.
## diet_premium = inclusion × (meal_price − displaced_cost)
## → displaced_cost = meal_price − diet_premium / inclusion
displaced_cost <- ulva_meal_price - diet_premium / ulva_inclusion

cat(sprintf("\n  Diet premium at current prices:    $%.2f/kg feed\n", diet_premium))
cat(sprintf("  Implied displaced ingredient cost: $%.2f/kg meal\n", displaced_cost))

### BAYESIAN MODELS

## MODEL 1: UNADJUSTED — diet + scaled start weight only 
# This is the optimistic ceiling and simulates potential benefit if feed rates were targeted by biomass
fit_unadjusted <- readRDS(here("models", "fit_weight_final_unadjusted.rds"))
summary(fit_unadjusted)

## MODEL 2: ADJUSTED — diet + per_capita_feed_z + start_ABW_z (primary)
# Adjusting for per-capita feed removes the feed-availability confound (Ulva tanks received more feed)
# Adjusting for start_ABW removes the baseline-size imbalance 
fit_adjusted <- readRDS(here("models", "fit_weight_final.rds"))
summary(fit_adjusted)

### EXTRACT POSTERIORS

b_unadjusted <- fit_unadjusted |>
  as_draws_df() |>
  pull(b_dietulva)

b_adjusted <- fit_adjusted |>
  as_draws_df() |>
  pull(b_dietulva)

growth_unadjusted <- exp(b_unadjusted) - 1
growth_adjusted   <- exp(b_adjusted)   - 1

cat(sprintf("\n  Growth — unadjusted (optimistic):  median %+.2f%%  P(>0) = %.3f\n",
            median(growth_unadjusted) * 100, mean(b_unadjusted > 0)))
cat(sprintf("  Growth — adjusted (primary):        median %+.2f%%  P(>0) = %.3f\n",
            median(growth_adjusted)   * 100, mean(b_adjusted   > 0)))

### PER-ANIMAL DATA FROM TRIAL 

## Get the control mean final weight (kg). This is the baseline animal whose growth is scaled in FW analysis

control_weight_kg <- df_econ |>
  filter(diet == "control") |>
  summarise(w = mean(weight_g) / 1000) |>
  pull(w)

## Per-capita feed in kg in control tank

feed_kg <- tank_df |>
  filter(diet == "control") |>
  summarise(f = mean(per_capita_feed) / 1e6) |>   # mg -> kg
  pull(f)

cat(sprintf("\n  Control mean final weight:         %.2f g\n", control_weight_kg * 1000))
cat(sprintf("  Per-animal feed (control basis):   %.2f g  (%.6f kg)\n",
            feed_kg * 1000, feed_kg))

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
days_ctrl <- log(TARGET_G / W0_ref) / ctrl_sgr     # days for control to reach target

cat(sprintf("\n  Common starting weight (control grand mean):  %.2f g\n", W0_ref))
cat(sprintf("  Control SGR:                                  %.4f %%/day\n", ctrl_sgr * 100))
cat(sprintf("  Days for control to reach %.0fg:              %.1f days\n", TARGET_G, days_ctrl))


### PART A — FINAL WEIGHT BREAK-EVEN (fixed time point, day 131)

## Question: Is the larger final weight in Ulva tanks enough to counteract diet premium?

## Derive economic functions

## Revenue gained per animal from Ulva relative to control:
revenue_gain <- function(growth) {
  control_weight_kg * growth * farmgate_price
}

## Extra feed cost per animal at a given Ulva meal price:
extra_feed_cost <- function(meal_price) {
  premium <- ulva_inclusion * (meal_price - displaced_cost)
  feed_kg * premium
}

## Net profit per animal = revenue gain − extra feed cost
net_profit <- function(growth, meal_price) {
  revenue_gain(growth) - extra_feed_cost(meal_price)
}

## Break-even Ulva meal price: the meal price at which net_profit = 0
breakeven_meal_price <- function(growth) {
  displaced_cost + revenue_gain(growth) / (feed_kg * ulva_inclusion)
}

## Summary for Part A:

summarise_economics <- function(growth, label) {
  be  <- breakeven_meal_price(growth)
  net <- net_profit(growth, ulva_meal_price)
  tibble(
    model             = label,
    growth_median_pct = round(median(growth) * 100, 2),
    rev_gain_mean     = round(mean(revenue_gain(growth)), 3),
    net_profit_mean   = round(mean(net), 3),
    net_profit_lo95   = round(quantile(net, 0.025), 3),
    net_profit_hi95   = round(quantile(net, 0.975), 3),
    p_profitable      = round(mean(net > 0), 3),
    be_meal_median    = round(median(be), 2),
    be_meal_lo95      = round(quantile(be, 0.025), 2),
    be_meal_hi95      = round(quantile(be, 0.975), 2)
  )
}

econ_summary <- bind_rows(
  summarise_economics(growth_unadjusted, "unadjusted"),
  summarise_economics(growth_adjusted,   "adjusted")
)

cat("\nPART A — FINAL WEIGHT ECONOMIC SUMMARY AT ACTUAL ULVA MEAL PRICE ($8.22/kg):\n")
print(as.data.frame(econ_summary), row.names = FALSE)

## Probability of profit curve

meal_grid_fw <- seq(0, 60, by = 0.5)

prob_curve <- function(growth, label) {
  tibble(meal_price = meal_grid_fw) |>
    rowwise() |>
    mutate(
      p_profitable = mean(net_profit(growth, meal_price) > 0),
      model        = label
    ) |>
    ungroup()
}

curve_df <- bind_rows(
  prob_curve(growth_unadjusted, "unadjusted"),
  prob_curve(growth_adjusted,   "adjusted")
)

p_breakeven <- ggplot(curve_df,
                      aes(x = meal_price, y = p_profitable,
                          colour = model, linetype = model)) +
  geom_vline(xintercept = ulva_meal_price,
             colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0.5,
             colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 1) +
  annotate("text",
           x = ulva_meal_price + 1, y = 0.05,
           label = sprintf("Actual meal\nprice $%.2f/kg", ulva_meal_price),
           hjust = 0, size = 3, colour = "#993C1D") +
  scale_colour_manual(
    values = model_cols,
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight")
  ) +
  scale_linetype_manual(
    values = c("unadjusted" = "solid", "adjusted" = "dashed"),
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight")
  ) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_x_continuous(labels = dollar_format()) +
  labs(
    x        = "Ulva meal price ($/kg)",
    y        = "P(Ulva profitable)"
  ) +
  theme_ulva() +
  theme(
    legend.position     = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background   = element_rect(fill = "white", colour = "grey85",
                                       linewidth = 0.3),
    legend.margin       = margin(4, 6, 4, 6)
  )

print(p_breakeven)
ggsave(here("figures", "p_breakeven.png"),
       plot = p_breakeven, dpi = 300, width = 9, height = 6, units = "in")


### PART B: DIET PREMIUM AS A FUNCTION OF ULVA MEAL PRICE

## Question: How does diet premium scale with price of raw Ulva meal?

premium_from_meal <- function(meal_price) {
  ulva_inclusion * (meal_price - displaced_cost)
}

premium_per_dollar <- ulva_inclusion

cat(sprintf("\n  Each $1/kg rise in Ulva meal price adds $%.2f/kg to the diet premium\n",
            premium_per_dollar))
cat(sprintf("  At actual meal price ($%.2f/kg): diet premium = $%.2f/kg feed\n",
            ulva_meal_price, premium_from_meal(ulva_meal_price)))

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
ggsave(here("figures", "p_premium.png"),
       plot = p_premium, dpi = 300, width = 8, height = 5, units = "in")

### PART C: TIME-TO-HARVEST BREAK-EVEN USING SGR

## Functions

# Ulva final weight per posterior draw, projected from W0_ref
ulva_Wf <- function(b) ctrl_Wf * exp(b)

# Ulva SGR per draw
ulva_sgr <- function(b) (log(ulva_Wf(b)) - log(W0_ref)) / DAYS

# Days for Ulva to reach target per draw
days_ulva <- function(b) log(TARGET_G / W0_ref) / ulva_sgr(b)

# Days saved relative to control per draw
days_saved <- function(b) days_ctrl - days_ulva(b)

# Cumulative feed to reach target (kg) for a given SGR
feed_to_target_kg <- function(sgr_val) {
  FEED_RATE * (TARGET_G - W0_ref) / 1000 / sgr_val
}

# Ulva diet price as a linear function of Ulva meal price (fixed inclusion)
ulva_diet_from_meal <- function(meal_price) {
  ctrl_diet_price + ulva_inclusion * (meal_price - displaced_cost)
}

# Feed cost to reach target for a given SGR and diet price
feed_cost <- function(sgr_val, diet_price) {
  feed_to_target_kg(sgr_val) * diet_price
}

# Net feed cost saving = control cost - Ulva cost (negative = Ulva costs more)
feed_saving <- function(b, meal_price) {
  dp <- ulva_diet_from_meal(meal_price)
  feed_cost(ctrl_sgr, ctrl_diet_price) - feed_cost(ulva_sgr(b), dp)
}

# Break-even Ulva diet price: diet price at which feed saving = 0
be_diet_price <- function(b) {
  feed_cost(ctrl_sgr, ctrl_diet_price) / feed_to_target_kg(ulva_sgr(b))
}

# Break-even Ulva meal price: backed out from the diet price break-even
be_meal_price <- function(b) {
  displaced_cost + (be_diet_price(b) - ctrl_diet_price) / ulva_inclusion
}

# Summary table
feed_cost_ctrl <- feed_cost(ctrl_sgr, ctrl_diet_price)

summarise_sgr_economics <- function(b, label) {
  tibble(
    model                  = label,
    sgr_adv_median_pct_day = round(median(b / DAYS * 100), 5),
    days_to_target_ctrl    = round(days_ctrl, 1),
    days_to_target_ulva    = round(median(days_ulva(b)), 1),
    days_saved_median      = round(median(days_saved(b)), 1),
    days_saved_lo95        = round(quantile(days_saved(b), 0.025), 1),
    days_saved_hi95        = round(quantile(days_saved(b), 0.975), 1),
    p_days_saved_gt0       = round(mean(days_saved(b) > 0), 3),
    feed_cost_ctrl         = round(feed_cost_ctrl, 4),
    feed_cost_ulva_median  = round(median(feed_cost(ulva_sgr(b), ulva_diet_price)), 4),
    net_saving_median      = round(median(feed_saving(b, ulva_meal_price)), 4),
    net_saving_lo95        = round(quantile(feed_saving(b, ulva_meal_price), 0.025), 4),
    net_saving_hi95        = round(quantile(feed_saving(b, ulva_meal_price), 0.975), 4),
    p_saving_gt0           = round(mean(feed_saving(b, ulva_meal_price) > 0), 3),
    be_meal_median         = round(median(be_meal_price(b)), 2),
    be_meal_lo95           = round(quantile(be_meal_price(b), 0.025), 2),
    be_meal_hi95           = round(quantile(be_meal_price(b), 0.975), 2),
    be_diet_median         = round(median(be_diet_price(b)), 3),
    be_diet_lo95           = round(quantile(be_diet_price(b), 0.025), 3),
    be_diet_hi95           = round(quantile(be_diet_price(b), 0.975), 3)
  )
}

sgr_summary <- bind_rows(
  summarise_sgr_economics(b_unadjusted, "unadjusted"),
  summarise_sgr_economics(b_adjusted,   "adjusted")
)

# PART C — SGR-BASED ECONOMIC SUMMARY AT ACTUAL ULVA MEAL PRICE ($8.22/kg):\n")
print(as.data.frame(sgr_summary), row.names = FALSE)

cat(sprintf("\n  Control feed cost to %.0fg:   $%.4f/animal\n", TARGET_G, feed_cost_ctrl))
cat(sprintf("  Current Ulva meal price:      $%.2f/kg\n\n", ulva_meal_price))

for (lbl in c("unadjusted", "adjusted")) {
  row <- sgr_summary |> filter(model == lbl)
  cat(sprintf("  [%s]\n", lbl))
  cat(sprintf("    SGR advantage:     median %+.5f %%/day\n", row$sgr_adv_median_pct_day))
  cat(sprintf("    Days saved:        median %.1f  (95%% CrI %.1f to %.1f)  P(>0) = %.3f\n",
              row$days_saved_median, row$days_saved_lo95, row$days_saved_hi95,
              row$p_days_saved_gt0))
  cat(sprintf("    Net feed saving:   median $%+.4f  (95%% CrI $%+.4f to $%+.4f)  P(>0) = %.3f\n",
              row$net_saving_median, row$net_saving_lo95, row$net_saving_hi95,
              row$p_saving_gt0))
  cat(sprintf("    Break-even meal:   median $%.2f/kg  (95%% CrI $%.2f to $%.2f)\n",
              row$be_meal_median, row$be_meal_lo95, row$be_meal_hi95))
  cat(sprintf("    Break-even diet:   median $%.3f/kg  (95%% CrI $%.3f to $%.3f)\n\n",
              row$be_diet_median, row$be_diet_lo95, row$be_diet_hi95))
}

## PART C: PROBABILITY-OF-POSITIVE-SAVING CURVE

# Create grid with different Ulva prices
meal_grid_sgr <- seq(0, 20, by = 0.25)

prob_saving_curve <- function(b, label) {
  tibble(meal_price = meal_grid_sgr) |>
    rowwise() |>
    mutate(
      p_saving_pos  = mean(feed_saving(b, meal_price) > 0),
      median_saving = median(feed_saving(b, meal_price)),
      model         = label
    ) |>
    ungroup()
}

saving_curve_df <- bind_rows(
  prob_saving_curve(b_unadjusted, "unadjusted"),
  prob_saving_curve(b_adjusted,   "adjusted")
)

# Highest meal price where P still meets threshold
thresh <- saving_curve_df |>
  group_by(model) |>
  summarise(
    price_above_95 = {
      vals <- meal_price[p_saving_pos >= 0.95]
      if (length(vals) == 0) NA_real_ else max(vals)
    },
    price_above_50 = {
      vals <- meal_price[p_saving_pos >= 0.50]
      if (length(vals) == 0) NA_real_ else max(vals)
    },
    .groups = "drop"
  )

print(as.data.frame(thresh), row.names = FALSE)

p_saving_curve <- ggplot(saving_curve_df,
                         aes(x = meal_price, y = p_saving_pos,
                             colour = model, linetype = model)) +
  geom_vline(xintercept = ulva_meal_price,
             colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0.95,
             colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_hline(yintercept = 0.50,
             colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 1) +
  annotate("text",
           x      = ulva_meal_price,
           y      = 0.12,
           label  = sprintf("Actual meal\nprice $%.2f/kg", ulva_meal_price),
           hjust  = -1,
           vjust  = 0,
           size   = 3,
           colour = "#993C1D") +
  scale_colour_manual(
    values = model_cols,
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight")
  ) +
  scale_linetype_manual(
    values = c("unadjusted" = "solid", "adjusted" = "dashed"),
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight")
  ) +
  scale_x_continuous(labels = dollar_format(), breaks = seq(0, 20, 2)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    x     = "Ulva meal price ($/kg)",
    y     = "P(feed cost saving > 0)",
    title = ""
  ) +
  theme_ulva() +
  theme(
    legend.position     = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background   = element_rect(fill = "white", colour = "grey85",
                                       linewidth = 0.3),
    legend.margin       = margin(4, 6, 4, 6)
  )

print(p_saving_curve)
ggsave(here("figures", "p_sgr_saving_by_meal_price.png"),
       plot = p_saving_curve, dpi = 300, width = 9, height = 6, units = "in")

### PART C: POSTERIOR DISTRIBUTION OF DAYS SAVED

days_df <- bind_rows(
  tibble(days_saved = days_saved(b_unadjusted), model = "unadjusted"),
  tibble(days_saved = days_saved(b_adjusted),   model = "adjusted")
)

p_days <- ggplot(days_df, aes(x = days_saved, fill = model, colour = model)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins      = 60,
                 alpha     = 0.45,
                 position  = "identity",
                 linewidth = 0.2) +
  geom_density(alpha = 0, linewidth = 0.8) +
  geom_vline(
    data = days_df |>
      group_by(model) |>
      summarise(med = median(days_saved), .groups = "drop"),
    aes(xintercept = med, colour = model),
    linewidth = 0.7, linetype = "dashed"
  ) +
  scale_fill_manual(
    values = model_cols,
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight")
  ) +
  scale_colour_manual(
    values = model_cols,
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight")
  ) +
  labs(
    x        = sprintf("Days saved to reach %.0fg harvest weight", TARGET_G),
    y        = "Density"
  ) +
  theme_ulva()

print(p_days)
ggsave(here("figures", "p_days_saved.png"),
       plot = p_days, dpi = 300, width = 9, height = 6, units = "in")

### EXPORT ALL DATA AS CSV

write.csv(econ_summary,    here("outputs", "ulva_economic_summary.csv"),        row.names = FALSE)
write.csv(curve_df,        here("outputs", "ulva_breakeven_curve.csv"),          row.names = FALSE)
write.csv(premium_df,      here("outputs", "ulva_premium_curve.csv"),            row.names = FALSE)
write.csv(sgr_summary,     here("outputs", "ulva_sgr_economic_summary.csv"),     row.names = FALSE)
write.csv(saving_curve_df, here("outputs", "ulva_sgr_saving_curve.csv"),         row.names = FALSE)

### END OF SCRIPT ###
