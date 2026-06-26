# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

## Step 2: Economic analysis

## Load required packages

## Run model using both the unadjusted (no feed corerection) and corrected effect size

# Set colour pallete
diet_cols <- c("unadjusted" = "#185FA5", "feed_adjusted" = "#0F6E56")

### ECONOMIC INPUTS 

## All prices in AUD - Feed quantities converted to kg.
farmgate_price   <- 25.00      # $/kg whole abalone at farmgate
ulva_diet_price  <- 3.78       # $/kg finished Ulva diet
ctrl_diet_price  <- 2.47       # $/kg finished control diet
ulva_meal_price  <- 8.22       # $/kg Ulva meal (the negotiable lever)
ulva_inclusion   <- 0.20       # Ulva meal inclusion level (20% )

diet_premium     <- ulva_diet_price - ctrl_diet_price   # $1.31/kg extra to feed Ulva

## The diet premium arises from including Ulva meal at 20%. Back out implied cost of ingredient(s) the meal displaces, so the premium 
## can be re-expressed as function of the Ulva meal e.g.   diet_premium = inclusion * (meal_price - displaced_cost)
displaced_cost   <- ulva_meal_price - diet_premium / ulva_inclusion

### PER-ANIMAL QUANTITIES FROM THE TRIAL DATA 

## Control mean final weight (kg) — the baseline animal whose growth is scaled.
## per_capita_feed is cumulative feed per animal over the 131-day trial, in mg;
## convert to kg. Use the CONTROL-tank feed as the cost basis: the higher Ulva
## feed reflects the (already removed) feed-availability confound, so charging it
## again on the cost side would double-count it.
control_weight_kg <- df |>
  filter(diet == "control") |>
  summarise(w = mean(weight_g) / 1000) |>
  pull(w)

feed_kg <- tank_df |>
  filter(diet == "control") |>
  summarise(f = mean(per_capita_feed) / 1e6) |>   # mg -> kg
  pull(f)

# Check
cat(sprintf("  Control mean final weight: %.2f g\n", control_weight_kg * 1000))
cat(sprintf("  Per-animal feed (control basis): %.2f g = %.6f kg\n",
            feed_kg * 1000, feed_kg))


## POSTERIOR GROWTH EFFECTS 

## Pull the full posterior of the Ulva effect from each weight model and convert
## from the log scale to a proportional weight effect via exp(beta) - 1. These
## vectors of draws — not point estimates — are what propagate through the
## economics, so the break-even comes out as a distribution.
growth_unadjusted <- fit_weight_informative_uncorrected |>   # diet + start_ABW_z 
  as_draws_df() |>
  mutate(growth = exp(b_dietulva) - 1) |>
  pull(growth)

growth_feed_adjusted <- fit_weight_informative |>      # diet + start_ABW_z + per_capita_feed_z
  as_draws_df() |>
  mutate(growth = exp(b_dietulva) - 1) |>
  pull(growth)

cat(sprintf("\n  Growth (unadjusted):   mean %.2f%%  P(>0) = %.3f\n",
            mean(growth_unadjusted) * 100, mean(growth_unadjusted > 0)))
cat(sprintf("  Growth (feed-adjusted): mean %.2f%%  P(>0) = %.3f\n",
            mean(growth_feed_adjusted) * 100, mean(growth_feed_adjusted > 0)))

### BREAK-EVEN FUNCTIONS 

## Per animal, over the trial:
##   revenue gain = control_weight * growth * farmgate_price   (a posterior vector)
##   extra feed cost at a given meal price m:
##       premium(m) = inclusion * (m - displaced_cost)
##       cost(m)    = feed_kg * premium(m)
##   net profit(m)  = revenue gain - cost(m)
##
## The break-even meal price for each posterior draw is the m where net = 0:
##   m_break = displaced_cost + (revenue gain) / (feed_kg * inclusion)

revenue_gain <- function(growth) control_weight_kg * growth * farmgate_price

breakeven_meal_price <- function(growth) {
  displaced_cost + revenue_gain(growth) / (feed_kg * ulva_inclusion)
}

net_profit <- function(growth, meal_price) {
  premium <- ulva_inclusion * (meal_price - displaced_cost)
  revenue_gain(growth) - feed_kg * premium
}


### SUMMARY AT THE ACTUAL PRICE 

summarise_economics <- function(growth, label) {
  be  <- breakeven_meal_price(growth)
  net <- net_profit(growth, ulva_meal_price)
  tibble(
    model            = label,
    rev_gain_mean    = mean(revenue_gain(growth)),
    net_profit_mean  = mean(net),
    net_profit_lo95  = quantile(net, 0.025),
    net_profit_hi95  = quantile(net, 0.975),
    p_profitable     = mean(net > 0),
    be_meal_median   = median(be),
    be_meal_lo95     = quantile(be, 0.025),
    be_meal_hi95     = quantile(be, 0.975)
  ) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

econ_summary <- bind_rows(
  summarise_economics(growth_unadjusted,    "unadjusted"),
  summarise_economics(growth_feed_adjusted, "feed_adjusted")
)

# Print summary
print(as.data.frame(econ_summary), row.names = FALSE)


### PROBABILITY-OF-PROFIT CURVE 
## For a grid of candidate Ulva meal prices, compute the posterior probability
## that Ulva is profitable (net > 0) under each growth model. This is the
## decision-support deliverable.
meal_grid <- seq(0, 60, by = 0.5)

prob_curve <- function(growth, label) {
  tibble(meal_price = meal_grid) |>
    rowwise() |>
    mutate(
      p_profitable = mean(net_profit(growth, meal_price) > 0),
      model = label
    ) |>
    ungroup()
}

curve_df <- bind_rows(
  prob_curve(growth_unadjusted,    "unadjusted"),
  prob_curve(growth_feed_adjusted, "feed_adjusted")
)


## PLOT 

p_breakeven <- ggplot(curve_df, aes(x = meal_price, y = p_profitable,
                                     colour = model, linetype = model)) +
  geom_vline(xintercept = ulva_meal_price, colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 1) +
  annotate("text", x = ulva_meal_price + 1, y = 0.05,
           label = "Actual meal price\n$8.22/kg", hjust = 0, size = 3, colour = "#993C1D") +
  scale_colour_manual(values = diet_cols, name = "Growth model") +
  scale_linetype_manual(values = c("unadjusted" = "solid",
                                    "feed_adjusted" = "dashed"), name = "Growth model") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_x_continuous(labels = dollar_format()) +
  labs(
    x = "Ulva meal price ($/kg)",
    y = "P(Ulva profitable)",
    title = "Probability Ulva is profitable vs Ulva meal price",
    subtitle = "Per-animal break-even, mortality excluded; full posterior propagated"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

print(p_breakeven)

# ggsave("ulva_breakeven_curve.png", p_breakeven, width = 8, height = 5, dpi = 300)


### EXPORT #####################################################################

setwd("C:/Users/RebeccaPedler/OneDrive - Yumbah")
write.csv(econ_summary, "ulva_economic_summary.csv", row.names = FALSE)
write.csv(curve_df,     "ulva_breakeven_curve.csv",  row.names = FALSE)

### DIET PREMIUM AS A FUNCTION OF ULVA MEAL PRICE #############################

## The diet premium is what actually hits the cost side; the meal price is the
## lever. Since Ulva is included at a fixed 20%, the premium scales linearly with
## the meal price: each $1/kg rise in the meal price adds (inclusion x $1) to the
## diet premium. The slope is therefore constant at $0.20/kg of premium per $1/kg
## of meal.

premium_from_meal <- function(meal_price) {
  ulva_inclusion * (meal_price - displaced_cost)
}

premium_per_dollar <- ulva_inclusion   # slope: $ premium added per $1 meal rise

cat(sprintf("\n  Each $1/kg rise in Ulva meal price adds $%.2f/kg to the diet premium\n",
            premium_per_dollar))
cat(sprintf("  At the actual meal price ($%.2f/kg), the diet premium is $%.2f/kg\n",
            ulva_meal_price, premium_from_meal(ulva_meal_price)))

## Premium across a range of candidate meal prices
premium_df <- tibble(meal_price = meal_grid) |>
  mutate(
    diet_premium   = premium_from_meal(meal_price),
    ulva_diet_cost = ctrl_diet_price + diet_premium   # implied finished-diet price
  )

# Show a readable subset
premium_df |>
  filter(meal_price %in% seq(0, 60, by = 10)) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  as.data.frame() |>
  print(row.names = FALSE)

## Plot: diet premium vs Ulva meal price, with the actual price marked
p_premium <- ggplot(premium_df, aes(x = meal_price, y = diet_premium)) +
  geom_vline(xintercept = ulva_meal_price, colour = "#993C1D", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 1, colour = "#0F6E56") +
  annotate("text", x = ulva_meal_price + 1, y = min(premium_df$diet_premium),
           label = sprintf("Actual meal price\n$%.2f/kg", ulva_meal_price),
           hjust = 0, vjust = 0, size = 3, colour = "#993C1D") +
  scale_x_continuous(labels = dollar_format()) +
  scale_y_continuous(labels = dollar_format()) +
  labs(
    x = "Ulva meal price ($/kg)",
    y = "Diet premium over control ($/kg feed)",
    title = "Diet premium added with each $1 rise in Ulva meal price",
    subtitle = sprintf("Slope fixed by %.0f%% inclusion: $%.2f premium per $1/kg meal",
                       ulva_inclusion * 100, premium_per_dollar)
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

print(p_premium)
