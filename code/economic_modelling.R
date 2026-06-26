# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 2: Economic analysis

## Load required packages

library(tidyverse)
library(brms)
library(here)
library(scales)  

# Set colour palette

model_cols <- c(
  "unadjusted" = "#185FA5",
  "adjusted"   = "#0F6E56"
)

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


### PRIORS

priors_logwt <- c(
  prior(student_t(3, 3.7, 1), class = Intercept),
  prior(normal(0, 1),         class = b),
  prior(exponential(1),       class = sigma)
)

### MODEL 1: UNADJUSTED — diet only (optimistic ceiling)

fit_unadjusted <- brm(
  formula = mean_log_weight ~ diet,
  data    = tank_df,
  family  = gaussian(),
  prior   = priors_logwt,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 42,
  control = list(adapt_delta = 0.95)
)

summary(fit_unadjusted)
pp_check(fit_unadjusted)

### MODEL 2: ADJUSTED — diet + per_capita_feed_z + start_ABW_z (primary)
fit_adjusted <- brm(
  formula = mean_log_weight ~ diet + per_capita_feed_z + start_ABW_z,
  data    = tank_df,
  family  = gaussian(),
  prior   = priors_logwt,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 42,
  control = list(adapt_delta = 0.95)
)

summary(fit_adjusted)
pp_check(fit_adjusted)

### POSTERIOR GROWTH EFFECTS

growth_unadjusted <- fit_unadjusted |>
  as_draws_df() |>
  mutate(growth = exp(b_dietulva) - 1) |>
  pull(growth)

growth_adjusted <- fit_adjusted |>
  as_draws_df() |>
  mutate(growth = exp(b_dietulva) - 1) |>
  pull(growth)

cat(sprintf("\n  Growth — unadjusted (optimistic): %+.2f%%  P(>0) = %.3f\n",
            median(growth_unadjusted) * 100, mean(growth_unadjusted > 0)))
cat(sprintf("  Growth — adjusted (primary):        %+.2f%%  P(>0) = %.3f\n",
            median(growth_adjusted)   * 100, mean(growth_adjusted   > 0)))

### PER-ABALONE PARAMETERS FROM TRIAL DATA

## Control mean final weight (kg) (using the arithmetic mean of individual-level weights)

control_weight_kg <- df_econ |>
  filter(diet == "control") |>
  summarise(w = mean(weight_g) / 1000) |>
  pull(w)

## Per-capita feed in kg using control-tank mean as cost basis (adjusted model already removes feed-availability confound)

feed_kg <- tank_df |>
  filter(diet == "control") |>
  summarise(f = mean(per_capita_feed) / 1e6) |>  # mg -> kg
  pull(f)

cat(sprintf("\n  Control mean final weight:         %.2f g\n", control_weight_kg * 1000))
cat(sprintf("  Per-animal feed (control tanks):   %.2f g  (%.6f kg)\n",
            feed_kg * 1000, feed_kg))

### ECONOMIC INPUTS

farmgate_prices <- c(20, 25, 30)   # $/kg whole abalone at several tangeable farmgates (considering ABARES and AAGA published values)
ulva_diet_price <-  3.78           # $/kg finished Ulva diet
ctrl_diet_price <-  2.47           # $/kg finished control diet
ulva_meal_price <-  8.22           # $/kg Ulva meal (excl. GST and delivered)
ulva_inclusion  <-  0.20           # Ulva meal inclusion rate (20%)

## Diet premium (extra cost of UL20 relative to control)
diet_premium <- ulva_diet_price - ctrl_diet_price  
diet_premium

## Back out cost of ingredient Ulva meal displaces.
## diet_premium = inclusion × (meal_price − displaced_cost) → displaced_cost = meal_price − diet_premium / inclusion
displaced_cost <- ulva_meal_price - diet_premium / ulva_inclusion

cat(sprintf("\n  Diet premium at current prices:    $%.2f/kg feed\n", diet_premium))
cat(sprintf("  Implied displaced ingredient cost: $%.2f/kg meal\n", displaced_cost))

### ECONOMIC FUNCTIONS

## Revenue gained per animal from Ulva relative to control:
##   control_weight × growth_proportion × farmgate_price
revenue_gain <- function(growth, farmgate_price) {
  control_weight_kg * growth * farmgate_price
}

## Extra feed cost per animal at a given Ulva meal price m:
##   premium(m) = inclusion × (m − displaced_cost)
##   cost(m)    = feed_kg × premium(m)

extra_feed_cost <- function(meal_price) {
  premium <- ulva_inclusion * (meal_price - displaced_cost)
  feed_kg * premium
}

## Net profit per animal = revenue gain − extra feed cost

net_profit <- function(growth, meal_price, farmgate_price) {
  revenue_gain(growth, farmgate_price) - extra_feed_cost(meal_price)
}

## Break-even Ulva meal price: the m at which net_profit = 0
##   m_break = displaced_cost + revenue_gain / (feed_kg × inclusion)

breakeven_meal_price <- function(growth, farmgate_price) {
  displaced_cost + revenue_gain(growth, farmgate_price) / (feed_kg * ulva_inclusion)
}

### SUMMARY AT THE ACTUAL ULVA MEAL PRICE

summarise_economics <- function(growth, label, farmgate_price) {
  be  <- breakeven_meal_price(growth, farmgate_price)
  net <- net_profit(growth, ulva_meal_price, farmgate_price)
  tibble(
    model             = label,
    farmgate_price    = farmgate_price,
    growth_median_pct = round(median(growth) * 100, 2),
    rev_gain_mean     = round(mean(revenue_gain(growth, farmgate_price)), 3),
    net_profit_mean   = round(mean(net), 3),
    net_profit_lo95   = round(quantile(net, 0.025), 3),
    net_profit_hi95   = round(quantile(net, 0.975), 3),
    p_profitable      = round(mean(net > 0), 3),
    be_meal_median    = round(median(be), 2),
    be_meal_lo95      = round(quantile(be, 0.025), 2),
    be_meal_hi95      = round(quantile(be, 0.975), 2)
  )
}

### BREAK-EVEN ANALYSIS BUNDLED BY FARMGATE PRICE

## Runs the full per-animal break-even analysis, plus minimum growth required to break-even with 99% certainty, at each farmgate price:

run_breakeven <- function(farmgate_price) {

  cat(sprintf("\n================  Farmgate price $%.2f/kg  ================\n",
              farmgate_price))

  ## Summary table for both growth models
  econ_summary <- bind_rows(
    summarise_economics(growth_unadjusted, "unadjusted", farmgate_price),
    summarise_economics(growth_adjusted,   "adjusted",   farmgate_price)
  )
  print(as.data.frame(econ_summary), row.names = FALSE)

  ## Minimum growth required at the current Ulva meal price for certain profitability
  cost_fixed  <- feed_kg * ulva_inclusion * (ulva_meal_price - displaced_cost)
  g_min       <- cost_fixed / (control_weight_kg * farmgate_price)

  p01_current <- quantile(growth_adjusted, 0.01) 
  shortfall   <- g_min - p01_current
  g_needed    <- median(growth_adjusted) + shortfall

  cat(sprintf("\n  Fixed extra feed cost per animal:        $%.4f\n", cost_fixed))
  cat(sprintf("  Break-even growth rate (certain):        %.4f (%.2f%%)\n",
              g_min, g_min * 100))
  cat(sprintf("  Growth needed for P(profitable) = 0.99:  %.4f (%.2f%%)\n",
              g_needed, g_needed * 100))
  cat(sprintf("  Current adjusted estimate (median):      %.4f (%.2f%%)\n",
              median(growth_adjusted), median(growth_adjusted) * 100))
  cat(sprintf("  Additional growth required:              %.4f (%.2f pp)\n",
              g_needed - median(growth_adjusted),
              (g_needed - median(growth_adjusted)) * 100))

  invisible(econ_summary)
}

### RUN ECONOMIC BREAK-EVEN ACROSS FARMGATE PRICE SCENARIOS

## Print the break-even analysis for each farmgate price

# Full table
econ_summary <- farmgate_prices |>
  map(run_breakeven) |>
  list_rbind()

# Summary
print(as.data.frame(econ_summary), row.names = FALSE)

### PROBABILITY-OF-PROFIT CURVE ACROSS FARMGATE PRICES

## Compute posterior probability across a range of Ulva meal prices, for each scenario (e.g., uncorrected vs corrected at each farmgate)

meal_grid <- seq(0, 60, by = 0.5)

prob_curve <- function(growth, label, farmgate_price) {
  tibble(meal_price = meal_grid) |>
    rowwise() |>
    mutate(
      p_profitable   = mean(net_profit(growth, meal_price, farmgate_price) > 0),
      model          = label,
      farmgate_price = farmgate_price
    ) |>
    ungroup()
}

curve_df <- farmgate_prices |>
  map(\(fg) bind_rows(
    prob_curve(growth_unadjusted, "unadjusted", fg),
    prob_curve(growth_adjusted,   "adjusted",   fg)
  )) |>
  list_rbind() |>
  mutate(farmgate_label = factor(sprintf("Farmgate $%.0f/kg", farmgate_price),
                                 levels = sprintf("Farmgate $%.0f/kg", farmgate_prices)))

# Plot break-even, one panel per farmgate price
p_breakeven <- ggplot(curve_df,
                      aes(x = meal_price, y = p_profitable,
                          colour = model, linetype = model)) +
  geom_vline(xintercept = ulva_meal_price,
             colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0.5,
             colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 1) +
  facet_wrap(~ farmgate_label) +
  scale_colour_manual(
    values = model_cols,
    name   = ""
  ) +
  scale_linetype_manual(
    values = c("unadjusted" = "solid", "adjusted" = "dashed"),
    name   = "Growth model",
    labels = c("unadjusted" = "Unadjusted (optimistic ceiling)",
               "adjusted"   = "Adjusted: feed + start weight (primary)")
  ) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_x_continuous(labels = dollar_format(accuracy = 0.01)) +
  labs(
    x        = expression(italic(Ulva) ~ "price ($/kg)"),
    y        = "P(Ulva profitable)",
    title    = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),  
    axis.ticks       = element_line(color = "black", size = 0.5), 
    legend.position  = "bottom"
  )

print(p_breakeven)
ggsave(here("figures", "p_breakeven_farmgate.png"), plot = p_breakeven, dpi = 300, width = 12, height = 5, units = "in")

### DIET PREMIUM AS A FUNCTION OF ULVA MEAL PRICE

premium_from_meal <- function(meal_price) {
  ulva_inclusion * (meal_price - displaced_cost)
}

premium_per_dollar <- ulva_inclusion

cat(sprintf("\n  Each $1/kg rise in Ulva meal price adds $%.2f/kg to the diet premium\n",
            premium_per_dollar))
cat(sprintf("  At actual meal price ($%.2f/kg): diet premium = $%.2f/kg feed\n",
            ulva_meal_price, premium_from_meal(ulva_meal_price)))

premium_df <- tibble(meal_price = meal_grid) |>
  mutate(
    diet_premium   = premium_from_meal(meal_price),
    ulva_diet_cost = ctrl_diet_price + diet_premium
  )

# Readable subset
cat("\nDiet premium at selected Ulva meal prices:\n")
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
  scale_x_continuous(labels = dollar_format(accuracy = 0.01)) +
  scale_y_continuous(labels = dollar_format(accuracy = 0.01)) +
  labs(
    x        = "Ulva meal price ($/kg)",
    y        = "Diet premium over control ($/kg feed)",
    title    = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),  
    axis.ticks       = element_line(color = "black", size = 0.5), 
    legend.position  = "bottom"
  )

# Print and save plot
print(p_premium)
ggsave(here("figures", "p_premium.png"), plot = p_premium, dpi = 300, width = 8, height = 5, units = "in")

### END OF SCRIPT ###
