# Heavy metals analysis

install.packges(c("tidyverse", "rstatix", "dunn.test"))

library(tidyverse)
library(rstatix)
library(dunn.test)

# Load data 

dat <- read.csv(here("data", "heavy_metals.csv"), stringsAsFactors = FALSE)

metals <- c("Arsenic", "Cadmium", "Lead", "Mercury", "Iron")

# Replace <LOD strings with 0.01, coerce to numeric
dat <- dat %>%
  mutate(across(all_of(metals),
                ~ as.numeric(ifelse(grepl("^<", trimws(.)), 0.01, .))))


## Create analysis function
# For a single metal × tissue combination:
#   1. Shapiro-Wilk normality (per group)
#   2. Levene homogeneity of variances
#   3. One-way ANOVA if both assumptions met; Kruskal-Wallis otherwise
#   4. Post-hoc: Tukey HSD (ANOVA) or Dunn with Bonferroni correction (KW)

run_analysis <- function(data, metal, tissue_type) {

  df <- data %>%
    filter(tissue == tissue_type) %>%
    select(diet, value = all_of(metal)) %>%
    filter(!is.na(value)) %>%
    mutate(diet = factor(diet))

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Metal:", metal, "| Tissue:", tissue_type, "\n")
  cat(strrep("-", 60), "\n")

  # Descriptive summary
  desc <- df %>%
    group_by(diet) %>%
    summarise(n    = n(),
              mean = mean(value),
              sd   = sd(value),
              min  = min(value),
              max  = max(value),
              .groups = "drop")
  print(desc)
  cat("\n")

  # Assumption 1: Shapiro-Wilk 
  sw <- df %>%
    group_by(diet) %>%
    summarise(
      shapiro_W = tryCatch(shapiro.test(value)$statistic, error = function(e) NA),
      shapiro_p = tryCatch(shapiro.test(value)$p.value,   error = function(e) NA),
      .groups = "drop"
    )
  cat("Shapiro-Wilk normality test:\n")
  print(sw)
  normality_ok <- all(sw$shapiro_p > 0.05, na.rm = TRUE)
  cat("  -> Normality assumption", ifelse(normality_ok, "MET", "VIOLATED"), "\n\n")

  # Assumption 2: Levene 
  lev <- levene_test(df, value ~ diet)
  cat("Levene test for homogeneity of variances:\n")
  print(lev)
  variance_ok <- lev$p > 0.05
  cat("  -> Equal-variances assumption", ifelse(variance_ok, "MET", "VIOLATED"), "\n\n")

  # Parametric vs non-parametric decision 
  use_parametric <- normality_ok && variance_ok

  if (use_parametric) {

    cat("Decision: ANOVA (both assumptions met)\n\n")

    aov_fit <- aov(value ~ diet, data = df)
    cat("One-way ANOVA:\n")
    print(summary(aov_fit))

    aov_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]

    if (!is.na(aov_p) && aov_p < 0.05) {
      cat("\nTukey HSD post-hoc:\n")
      print(TukeyHSD(aov_fit))
    } else {
      cat("\nANOVA non-significant (p =", round(aov_p, 4),
          "); post-hoc not performed.\n")
    }

    test_stat <- summary(aov_fit)[[1]][["F value"]][1]
    p_value   <- aov_p

  } else {

    cat("Decision: Kruskal-Wallis (assumption(s) violated)\n\n")

    kw <- kruskal.test(value ~ diet, data = df)
    cat("Kruskal-Wallis test:\n")
    print(kw)

    if (kw$p.value < 0.05) {
      cat("\nDunn post-hoc (Bonferroni correction):\n")
      dunn.test(df$value, df$diet, method = "bonferroni", kw = FALSE, label = TRUE)
    } else {
      cat("\nKruskal-Wallis non-significant (p =",
          round(kw$p.value, 4), "); post-hoc not performed.\n")
    }

    test_stat <- kw$statistic
    p_value   <- kw$p.value
  }

  # Return summary row
  data.frame(
    Metal        = metal,
    Tissue       = tissue_type,
    Test         = ifelse(use_parametric, "ANOVA", "Kruskal-Wallis"),
    Normality    = ifelse(normality_ok, "Met", "Violated"),
    Eq_Variances = ifelse(variance_ok,  "Met", "Violated"),
    Test_stat    = round(unname(test_stat), 3),
    p_value      = round(p_value, 4),
    Significant  = ifelse(p_value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}


# Run across all metals × tissues 

tissues <- c("viscera", "foot")
results <- list()

for (tis in tissues) {
  for (met in metals) {
    results[[paste(met, tis, sep = "_")]] <- run_analysis(dat, met, tis)
  }
}

# Summary table
 
summary_table <- do.call(rbind, results)
rownames(summary_table) <- NULL
 
cat("\n", strrep("=", 60), "\n", sep = "")
cat("SUMMARY TABLE\n")
cat(strrep("=", 60), "\n\n")
print(summary_table, row.names = FALSE)
