# Helper functions used across the modelling scripts

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
