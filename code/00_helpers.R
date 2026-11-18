# Project: A probabilistic cost–benefit analysis of macroalgal dietary supplementation in commercial greenlip abalone (Haliotis laevigata) aquaculture

# Part 00: Helper functions used across the modelling scripts

## FUNCTION FOR GETTING POSTERIOR PROBABILITIES
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

# FUNCTION FOR GETTING POSTERIOR PROBABILITY PLOTS

# Font
windowsFonts(`Times New Roman` = windowsFont("Times New Roman"))
                 
# Adapted from Recio et al. (2025), Animal Behaviour 222:123106
# (github.com/Pablo-Recio/CORT_Temp_Behavioural_flexibility, R/func.r:
#  df_plotBD1, df_plotBD2 and plotting(), which draw Figure 2b and 2d).

# Get posterior draws
get_effects <- function(model, specs, transform = identity, panel = NULL) {
  draws <- if (inherits(model, "brmsfit")) {
    as.data.frame(brms::as_draws_df(model), check.names = FALSE)
  } else {
    as.data.frame(model, check.names = FALSE)
  }
  exprs <- rlang::parse_exprs(unname(specs))
  vals  <- lapply(exprs, function(e) transform(rlang::eval_tidy(e, data = draws)))
  out <- tibble::tibble(
    term  = factor(rep(names(specs), each = nrow(draws)), levels = names(specs)),
    value = unlist(vals, use.names = FALSE)
  )
  if (!is.null(panel)) out$panel <- panel
  out
}
                  
# Numeric summary (same numbers as the plot; handy for a results table)
summarise_effects <- function(effects, ref = 0) {
  grp <- intersect(c("panel", "term"), names(effects))
  effects |>
    group_by(across(all_of(grp))) |>
    summarise(
      mean   = mean(value),
      median = median(value),
      sd     = sd(value),
      lower95 = quantile(value, 0.025),
      upper95 = quantile(value, 0.975),
      lower66 = quantile(value, 0.17),
      upper66 = quantile(value, 0.83),
      p_gt_ref = mean(value > ref),
      .groups = "drop"
    )
}

# Plot
plot_effects <- function(effects,
                         colours      = NULL,
                         labels       = NULL,
                         ref          = 0,
                         centre       = c("mean", "median"),
                         interval     = c("sd", "qi"),
                         xlab         = "Effect estimate",
                         slab_height  = 0.8,
                         slab_alpha   = 0.6,
                         show_prob    = FALSE,
                         base_size    = 12,
                         title_size   = 12,   
                         ylab_size    = 12,   
                         xtext_size   = 12,   
                         facet_scales = "free_x",
                         facet_ncol   = NULL) {

  centre   <- match.arg(centre)
  interval <- match.arg(interval)
  has_panel <- "panel" %in% names(effects)
  grp   <- if (has_panel) c("panel", "term") else "term"
  terms <- levels(effects$term)
  n     <- length(terms)
  ypos  <- function(term) n - as.integer(term) + 1        

  dens <- effects |>
    group_by(across(all_of(grp))) |>
    group_modify(~ {
      d <- stats::density(.x$value, n = 512)
      tibble::tibble(x = d$x, dens = d$y / max(d$y))
    }) |>
    ungroup() |>
    mutate(y0 = ypos(term), y1 = y0 + slab_height * dens)

  summ <- summarise_effects(effects, ref) |>
    mutate(y0 = ypos(term),
           centre_val = if (centre == "mean") mean else median,
           lab = case_when(p_gt_ref > 0.995 ~ sprintf("P(> %g) > 0.99", ref),
                           p_gt_ref < 0.005 ~ sprintf("P(> %g) < 0.01", ref),
                           TRUE             ~ sprintf("P(> %g) = %.2f", ref, p_gt_ref)))
  if (interval == "sd") {
    summ <- mutate(summ, lo_thick = centre_val - sd, hi_thick = centre_val + sd)
  } else {
    summ <- mutate(summ, lo_thick = lower66, hi_thick = upper66)
  }

  if (is.null(colours)) colours <- stats::setNames(scales::hue_pal()(n), terms)
  lab_expr <- if (is.null(labels)) sprintf("bold('%s')", terms) else unname(labels[terms])

  p <- ggplot() +
    geom_vline(xintercept = ref, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_ribbon(data = dens,
                aes(x = x, ymin = y0, ymax = y1, fill = term, group = term),
                alpha = slab_alpha, colour = "grey20", linewidth = 0.3,
                outline.type = "upper")

  if (interval == "qi") {
    p <- p + geom_segment(data = summ,
                          aes(x = lower95, xend = upper95, y = y0, yend = y0),
                          linewidth = 0.6, colour = "black")
  }
  p <- p +
    geom_segment(data = summ,
                 aes(x = lo_thick, xend = hi_thick, y = y0, yend = y0),
                 linewidth = 1.8, colour = "black") +
    geom_point(data = summ, aes(x = centre_val, y = y0), size = 3, colour = "black")

  if (show_prob) {
    p <- p + geom_text(data = summ, aes(x = Inf, y = y0 + 0.12, label = lab), hjust = 1.05,
                       vjust = 0, size = base_size * 0.26, colour = "grey25", family = "Times New Roman")
  }

  p <- p +
    scale_fill_manual(values = colours, guide = "none") +
    scale_y_continuous(breaks = n:1, labels = parse(text = lab_expr),
                       expand = expansion(add = c(0.25, 0.15))) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid   = element_blank(),
      axis.line.x  = element_line(colour = "black", linewidth = 0.6),
      axis.line.y  = element_line(colour = "black", linewidth = 0.6),
      axis.text    = element_text(colour = "black"),
      axis.text.x  = element_text(size = xtext_size),
      axis.text.y  = element_text(size = ylab_size),
      axis.title.x = element_text(size = title_size, margin = margin(t = 10)),
      strip.text   = element_text(face = "bold"),
      text         = element_text(family = "Times New Roman")
    )
  if (has_panel) p <- p + facet_wrap(~ panel, scales = facet_scales, ncol = facet_ncol)
  p
}
