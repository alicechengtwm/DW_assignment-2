# IMPORTANT: the original dataset is not publicly distributable.
# This script generates a synthetic dataset according to the paper information as close as possible, 
# but it cannot reproduce the authors' exact estimates.

OUTPUT_DIR <- "case_study_output" #change to your output directory
N_IMPUTATIONS <- 45
MAX_ITERATIONS <- 20
RANDOM_SEED <- 12022025
SYNTHETIC_N <- 11506

library(mice)
library(nnet)

# Create variables
binary_vars <- c(
  "death1yr", "female", "acpulmed", "cshock", "diabetes", "highbp",
  "cva", "dyslip", "famhxcad", "angina", "cancer", "dementia", "pud",
  "prevmi", "asthma", "depres", "perartdis", "chf", "hyprthyr", "as",
  "prevrevasc"
)
continuous_vars <- c(
  "age", "sysbp", "diasbp", "hrtrate", "resp", "hgb", "wbc", "sod",
  "pot", "glucose", "urea", "cr"
)
smoking_levels <- c("NEVER", "CURRENT", "RECENT", "FORMER")
required_vars <- c(
  "death1yr", "female", "age", "acpulmed", "cshock", "diabetes", "highbp",
  "smokhx", "cva", "dyslip", "famhxcad", "angina", "cancer", "dementia",
  "pud", "prevmi", "asthma", "depres", "perartdis", "chf", "hyprthyr",
  "as", "prevrevasc", "sysbp", "diasbp", "hrtrate", "resp", "hgb", "wbc",
  "sod", "pot", "glucose", "urea", "cr"
)

# Create regression formula
analysis_formula <- death1yr ~ female + age + acpulmed + cshock +
  diabetes + highbp + cva + dyslip + famhxcad + angina + cancer +
  dementia + pud + prevmi + asthma + depres + perartdis + chf +
  hyprthyr + as + prevrevasc + sysbp + diasbp + hrtrate + resp +
  hgb + wbc + sod + pot + glucose + urea + cr + smokhx

# Create synthetic data function
make_synthetic_data <- function(n, seed) {
  set.seed(seed)

  # Generate the nine continuous variables highlighted in Section 3.2.1.
  # Means/SDs and correlations are clinically plausible assumptions because the
  # paper reports that these were estimated from the restricted data but does
  # not publish their numerical values.
  x_names <- c("age", "sysbp", "hrtrate", "resp", "glucose", "wbc", "urea", "cr", "hgb")
  x_means <- c(68, 135, 82, 20, 7.2, 8.0, 7.0, 100, 130)
  x_sds   <- c(13, 25, 18, 5, 3.0, 3.5, 4.0, 50, 20)
  cor_x <- matrix(c(
    1.00, -0.20, 0.10, 0.10, 0.10, 0.10, 0.25, 0.25, -0.20,
    -0.20, 1.00, -0.05, -0.05, -0.05, 0.00, -0.10, -0.10, 0.10,
    0.10, -0.05, 1.00, 0.25, 0.10, 0.20, 0.05, 0.05, -0.10,
    0.10, -0.05, 0.25, 1.00, 0.10, 0.20, 0.05, 0.05, -0.10,
    0.10, -0.05, 0.10, 0.10, 1.00, 0.25, 0.20, 0.10, -0.10,
    0.10, 0.00, 0.20, 0.20, 0.25, 1.00, 0.20, 0.15, -0.10,
    0.25, -0.10, 0.05, 0.05, 0.20, 0.20, 1.00, 0.45, -0.15,
    0.25, -0.10, 0.05, 0.05, 0.10, 0.15, 0.45, 1.00, -0.15,
    -0.20, 0.10, -0.10, -0.10, -0.10, -0.10, -0.15, -0.15, 1.00
  ), nrow = 9, byrow = TRUE)
  z <- matrix(rnorm(n * 9), nrow = n, ncol = 9) %*% chol(cor_x)
  x <- sweep(sweep(z, 2, x_sds, `*`), 2, x_means, `+`)
  colnames(x) <- x_names
  d <- as.data.frame(x)

  # Keep simulated measurements within plausible ranges.
  d$age <- pmin(pmax(d$age, 18), 100)
  d$sysbp <- pmin(pmax(d$sysbp, 60), 250)
  d$hrtrate <- pmin(pmax(d$hrtrate, 30), 180)
  d$resp <- pmin(pmax(d$resp, 5), 50)
  d$glucose <- pmin(pmax(d$glucose, 2), 35) # mmol/L
  d$wbc <- pmin(pmax(d$wbc, 1), 40)          # x10^9/L
  d$urea <- pmin(pmax(d$urea, 1), 50)        # mmol/L
  d$cr <- pmin(pmax(d$cr, 20), 900)          # umol/L
  d$hgb <- pmin(pmax(d$hgb, 50), 200)        # g/L
  d$diasbp <- pmin(pmax(rnorm(n, 78 - 0.10 * (d$age - 65), 14), 25), 150)
  d$sod <- pmin(pmax(rnorm(n, 138, 4), 110), 170)
  d$pot <- pmin(pmax(rnorm(n, 4.1, 0.5), 2), 8)

  # Age and sex are complete in the paper. Other binary prevalences are not
  # reported in the article; these are explicit synthetic assumptions.
  d$female <- rbinom(n, 1, plogis(qlogis(0.38) - 0.01 * (d$age - 68)))
  binary_prob <- c(
    acpulmed = 0.22, cshock = 0.06, diabetes = 0.24, highbp = 0.55,
    cva = 0.12, dyslip = 0.40, famhxcad = 0.28, angina = 0.30,
    cancer = 0.10, dementia = 0.04, pud = 0.08, prevmi = 0.25,
    asthma = 0.08, depres = 0.12, perartdis = 0.10, chf = 0.15,
    hyprthyr = 0.025, as = 0.06, prevrevasc = 0.18
  )
  for (v in names(binary_prob)) {
    age_slope <- if (v %in% c("cshock", "acpulmed")) 0.012 else 0.025
    d[[v]] <- rbinom(n, 1, plogis(qlogis(binary_prob[[v]]) + age_slope * (d$age - 68)))
  }

  # The paper reports observed smoking-history proportions of 33.8% never,
  # 38.0% current, 2.0% recent, and 26.3% former among records with observed
  # smoking history. Use these rounded proportions for the synthetic full data.
  d$smokhx <- sample(smoking_levels, n, replace = TRUE,
                     prob = c(0.338, 0.380, 0.020, 0.262))

  # Generate one-year mortality using reported Table 1 odds ratios for smoking
  # and selected comorbidities, plus explicitly assumed effects for variables
  # whose coefficients are not reproduced here. The intercept
  # is calibrated to the reported 20.0% mortality; this is not the paper's exact
  # fitted model because the original data and full coefficient vector are not public.
  lp_without_intercept <-
    log(1.071) * ((d$age - 65) / 1) +
    log(1.013) * d$acpulmed + log(5.141) * d$cshock +
    log(1.005) * d$diabetes + log(1.117) * d$highbp +
    log(1.437) * d$cva + log(0.839) * d$dyslip +
    log(0.920) * d$famhxcad + log(1.235) * d$angina +
    log(1.311) * d$cancer + log(1.605) * d$dementia +
    log(0.783) * d$pud + log(1.180) * d$prevmi +
    log(0.834) * d$asthma + log(1.323) * d$depres +
    log(1.267) * d$perartdis + log(1.458) * d$chf +
    log(0.880) * d$hyprthyr + log(1.041) * d$as +
    log(1.180) * d$prevrevasc +
    log(0.984) * ((d$sysbp - 135) / 10) +
    log(1.015) * ((d$diasbp - 78) / 10) +
    log(1.054) * ((d$hrtrate - 82) / 10) +
    log(1.041) * ((d$resp - 20) / 5) +
    log(1.041) * ((d$hgb - 130) / 10) +
    log(1.041) * ((d$wbc - 8) / 2) +
    log(1.028) * ((d$sod - 138) / 5) +
    log(1.054) * ((d$pot - 4.1) / 0.5) +
    log(1.041) * ((d$glucose - 7.2) / 2) +
    log(1.041) * ((d$urea - 7) / 3) +
    log(1.002) * (d$cr - 100) +
    ifelse(d$smokhx == "CURRENT", log(1.177),
    ifelse(d$smokhx == "FORMER", log(0.991),
    ifelse(d$smokhx == "RECENT", log(1.200), 0)))
  target_mortality <- 0.20
  intercept <- uniroot(function(a) mean(plogis(a + lp_without_intercept)) - target_mortality,
                       interval = c(-50, 20))$root
  d$death1yr <- rbinom(n, 1, plogis(intercept + lp_without_intercept))

  # Save the complete synthetic dataset before imposing missingness.
  write.csv(d, file.path(OUTPUT_DIR, "synthetic_complete_dataset.csv"), row.names = FALSE)

  # The paper reports 14.6% missing smoking history and 45.0% with at least
  # one missing variable. It does not publish the full variable-by-variable
  # missingness pattern. Enforce those two reported rates, then distribute
  # additional missing values over eligible predictors as a documented
  # synthetic assumption. Age, sex, and one-year mortality remain complete.
  n_any_missing <- round(0.45 * n)
  rows_any_missing <- sample.int(n, n_any_missing)
  n_smoking_missing <- round(0.146 * n)
  rows_smoking_missing <- sample(rows_any_missing, n_smoking_missing)
  d$smokhx[rows_smoking_missing] <- NA

  eligible_for_missing <- setdiff(required_vars, c("death1yr", "female", "age", "smokhx"))
  for (i in rows_any_missing) {
    # Each selected record receives at least one missing predictor. Additional
    # missingness produces a more realistic multivariable incomplete-data pattern.
    n_other_missing <- if (i %in% rows_smoking_missing) sample(0:2, 1) else sample(1:3, 1)
    if (n_other_missing > 0) {
      vars_i <- sample(eligible_for_missing, min(n_other_missing, length(eligible_for_missing)))
      d[i, vars_i] <- NA
    }
  }

  d$smokhx <- factor(d$smokhx, levels = smoking_levels)
  d
}

# Generate synthetic dataset
dat <- make_synthetic_data(SYNTHETIC_N, RANDOM_SEED)

names(dat) <- trimws(names(dat))
dat <- dat[, required_vars, drop = FALSE]
for (v in unique(c(binary_vars, continuous_vars))) {
  dat[[v]] <- suppressWarnings(as.numeric(as.character(dat[[v]])))
}
dat$smokhx <- factor(toupper(as.character(dat$smokhx)), levels = smoking_levels)
dat$smokhx <- relevel(dat$smokhx, ref = "NEVER")

write.csv(dat, file.path(OUTPUT_DIR, "synthetic_observed_dataset_with_missingness.csv"), row.names = FALSE)

# Check the missingness
missingness <- data.frame(
  variable = names(dat),
  missing_n = vapply(dat, function(x) sum(is.na(x)), numeric(1)),
  missing_percent = round(100 * vapply(dat, function(x) mean(is.na(x)), numeric(1)), 2)
)
write.csv(missingness, file.path(OUTPUT_DIR, "missingness_summary.csv"), row.names = FALSE)

