#### Preparations ----

# Read temperature data
temps <- read.delim("./data/era5_thermal_regimes.txt")

# Load functions
source("./R/functions.R")

#### Sensitivity 1: accuracy of the lambda = R0^(1/T) approximation ----

# Compare the lambda = R0^(1/T_gen) approximation of both models with exact
# Euler-Lotka solutions per cohort, for 50 randomized site x warming draws

# Solve sum(l(x) * m(x) * lambda^-x) = 1 for one cohort (log space for small lx)
euler_lotka_lambda <- function(x, mx, lx) {

  keep <- mx > 0
  if (length(keep) == 0 || lx <= 0 || !any(keep)) return(0)
  x <- x[keep]
  mx <- mx[keep]

  log_terms <- log(lx) + log(mx)
  f <- function(r) {
    s <- sum(exp(log_terms - r * x))
    if (!is.finite(s)) s <- 1e300
    s - 1
  }

  # Bracket around the R0^(1/T) approximation and expand until a sign change
  R0 <- sum(mx) * lx
  T_gen <- sum(x * mx) / sum(mx)
  r0 <- log(R0) / T_gen
  lo <- r0 - 0.1
  hi <- r0 + 0.1
  while (f(lo) < 0) lo <- lo - 0.1
  while (f(hi) > 0) hi <- hi + 0.1

  r <- uniroot(f, c(lo, hi), tol = 1e-12)$root
  return(exp(r))
}

# Mirrors composite_fitness_function() (functions.R), adding exact Euler-Lotka lambdas
composite_lambda_exact <- function(temperature, origin) {

  dev_r <- development_rate(temperature = temperature, origin = origin)
  growth_r <- growth_rate(temperature = temperature, origin = origin)
  surv_r <- survival_rate(temperature = temperature, origin = origin)
  fert_r <- fertility_rate(temperature = temperature, origin = origin)

  raw_dev_time <- 1 / dev_r
  T_dev_hours <- ceiling(raw_dev_time)
  T_dev_hours <- ifelse(!is.na(T_dev_hours) & T_dev_hours > 250 * 24, 250 * 24, T_dev_hours)

  mass_ind <- growth_r * T_dev_hours
  mass_ind <- ifelse(!is.na(mass_ind) & mass_ind < 0.0005, 0.0005, mass_ind)
  mass_ind <- ifelse(!is.na(mass_ind) & mass_ind > 0.0035, 0.0035, mass_ind)

  survival_total <- surv_r ^ T_dev_hours
  survival_total[is.na(survival_total)] <- 0
  fertility_total <- fert_r ^ T_dev_hours
  fertility_total[is.na(fertility_total)] <- 0

  max_fec <- mass_to_fecundity(mass = mass_ind, origin = origin) / 2

  R0 <- survival_total * fertility_total * max_fec

  ovi_fraction <- fecundity_to_oviposition_rate(1, temperature)
  raw_duration <- 1 / ovi_fraction
  k_full_hours <- floor(raw_duration)

  remainder_weight <- 1 - (k_full_hours * ovi_fraction)
  remainder_weight <- ifelse(remainder_weight < 0, 0, remainder_weight)

  moment_full <- ovi_fraction * (k_full_hours * (k_full_hours + 1) / 2)
  moment_remainder <- remainder_weight * (k_full_hours + 1)
  mean_ovi_time <- moment_full + moment_remainder
  T_gen_hours <- T_dev_hours + mean_ovi_time

  fitness_lambda_hourly <- R0^(1 / T_gen_hours)
  fitness_lambda_hourly[is.na(fitness_lambda_hourly)] <- 0

  # Exact solution per hourly cohort
  lambda_exact <- numeric(length(temperature))
  for (i in seq_along(temperature)) {
    x_hours <- T_dev_hours[i] + 1:(k_full_hours[i] + 1)
    mx <- max_fec[i] * c(rep(ovi_fraction[i], k_full_hours[i]), remainder_weight[i])
    lambda_exact[i] <- euler_lotka_lambda(x_hours, mx, survival_total[i] * fertility_total[i])
  }

  return(list(lambda.approx = fitness_lambda_hourly, lambda.exact = lambda_exact))
}

# Mirrors separate_fitness_function() (functions.R), adding exact Euler-Lotka lambdas
separate_lambda_exact <- function(temperatures, origin) {

  n_temps <- length(temperatures)
  main_cohorts <- 1:365
  main_start_hours <- (364 + main_cohorts) * 24 + 1
  if (max(main_start_hours) > n_temps) stop("Temperature vector too short.")

  all_dev_rates <- development_rate(temperatures, origin)
  all_growth_rates <- growth_rate(temperatures, origin)
  all_surv_rates <- survival_rate(temperatures, origin)
  all_fert_rates <- fertility_rate(temperatures, origin)

  C_lethal_count <- c(0, cumsum(all_surv_rates <= 0))
  C_infertile_count <- c(0, cumsum(all_fert_rates <= 0))

  all_surv_rates[all_surv_rates <= 0] <- 1e-16
  all_fert_rates[all_fert_rates <= 0] <- 1e-16

  all_surv_rates_log <- log(all_surv_rates)
  all_fert_rates_log <- log(all_fert_rates)
  ovi_temp_factors <- fecundity_to_oviposition_rate(1, temperatures)

  C_dev <- c(0, cumsum(all_dev_rates))
  C_growth <- c(0, cumsum(all_growth_rates))
  C_surv_log <- c(0, cumsum(all_surv_rates_log))
  C_fert_log <- c(0, cumsum(all_fert_rates_log))

  dev_targets <- 1 + C_dev[main_start_hours]
  natural_end_hours <- findInterval(dev_targets, C_dev)
  raw_duration <- natural_end_hours - main_start_hours + 1
  max_duration_hours <- 250 * 24
  simulated_dev_time_hours <- ifelse(is.na(raw_duration) | raw_duration > max_duration_hours,
                                     max_duration_hours,
                                     raw_duration)
  dev_end_hours <- main_start_hours + simulated_dev_time_hours - 1
  dev_end_hours[dev_end_hours >= n_temps] <- n_temps

  simulated_mass <- C_growth[dev_end_hours + 1] - C_growth[main_start_hours]
  simulated_mass[simulated_mass < 0.0005] <- 0.0005
  simulated_mass[simulated_mass > 0.0035] <- 0.0035

  total_survival <- exp(C_surv_log[dev_end_hours + 1] - C_surv_log[main_start_hours])
  lethal_hits <- C_lethal_count[dev_end_hours + 1] - C_lethal_count[main_start_hours]
  total_survival[lethal_hits > 0] <- 0

  total_fertility <- exp(C_fert_log[dev_end_hours + 1] - C_fert_log[main_start_hours])
  total_fertility[is.na(total_fertility)] <- 0
  infertile_hits <- C_infertile_count[dev_end_hours + 1] - C_infertile_count[main_start_hours]
  total_fertility[infertile_hits > 0] <- 0

  simulated_lambda <- numeric(365)
  lambda_exact <- numeric(365)

  for (i in 1:365) {

    total_egg_pool <- mass_to_fecundity(simulated_mass[i], origin) / 2

    adult_start_hour <- dev_end_hours[i] + 1
    eggs_laid_so_far <- 0
    hour_counter <- 0
    hourly_offspring_vector <- c()

    while (eggs_laid_so_far < total_egg_pool) {
      current_hour <- adult_start_hour + hour_counter
      if (current_hour > n_temps) break

      eggs_this_hour <- total_egg_pool * ovi_temp_factors[current_hour]
      # Clamp eggs to remaining pool
      eggs_remaining <- total_egg_pool - eggs_laid_so_far
      if(eggs_this_hour > eggs_remaining) eggs_this_hour <- eggs_remaining

      hourly_offspring_vector <- c(hourly_offspring_vector, eggs_this_hour)
      eggs_laid_so_far <- eggs_laid_so_far + eggs_this_hour
      hour_counter <- hour_counter + 1

      if(hour_counter > 200 * 24) break
    }

    age_of_mother_hours <- simulated_dev_time_hours[i] + (1:length(hourly_offspring_vector))
    mx <- hourly_offspring_vector
    x_hours <- age_of_mother_hours
    lx <- total_survival[i] * total_fertility[i]

    R0 <- sum(mx) * lx
    T_gen_hours <- sum(x_hours * mx) / sum(mx)
    simulated_lambda[i] <- R0 ^ (1 / (T_gen_hours))

    lambda_exact[i] <- euler_lotka_lambda(x_hours, mx, lx)
  }

  simulated_lambda[is.na(simulated_lambda)] <- 0

  return(list(lambda.approx = simulated_lambda, lambda.exact = lambda_exact))
}

##########

#### Run the 50 randomized site x warming draws ----

set.seed(1)

# Sites with at least three years of hourly data (cf. fitness_predictions.R)
site_counts <- table(temps$unique.ID)
eligible_sites <- names(site_counts)[site_counts >= 365 * 24 * 3 - 24]

sampled_sites <- sample(eligible_sites, 50)
sampled_warming <- sample(seq(from = 0, to = 5.5, by = 0.1), 50, replace = TRUE)

year_hours <- 365 * 24
results_list <- list()

for (draw in 1:50) {

  location_data <- subset(temps, unique.ID == sampled_sites[draw])
  location_name <- location_data$location[1]
  warmed_temperatures <- location_data$temperature + sampled_warming[draw]

  # Separate components model (daily cohorts of year 2)
  sep_res <- separate_lambda_exact(warmed_temperatures, location_name)

  # Composite model (hourly values of year 2, as in calculate_fitness)
  temps_year2 <- warmed_temperatures[(year_hours + 1):(2 * year_hours)]
  comp_res <- composite_lambda_exact(temps_year2, location_name)

  # Cross-check: the mirrored approximations must reproduce the pipeline output
  stopifnot(isTRUE(all.equal(sep_res$lambda.approx,
                             separate_fitness_function(warmed_temperatures, location_name)$lambda.hourly)))
  stopifnot(isTRUE(all.equal(comp_res$lambda.approx,
                             composite_fitness_function(temps_year2, location_name))))

  results_list[[length(results_list) + 1]] <- data.frame(
    unique.ID = sampled_sites[draw],
    location = location_name,
    warming = sampled_warming[draw],
    model = "separate",
    lambda.approx = sep_res$lambda.approx,
    lambda.exact = sep_res$lambda.exact
  )

  results_list[[length(results_list) + 1]] <- data.frame(
    unique.ID = sampled_sites[draw],
    location = location_name,
    warming = sampled_warming[draw],
    model = "composite",
    lambda.approx = comp_res$lambda.approx,
    lambda.exact = comp_res$lambda.exact
  )

  cat(paste("  - Completed draw", draw, "of 50:", sampled_sites[draw], "+", sampled_warming[draw], "C\n"))
}

lambda_check <- do.call(rbind, results_list)

# Summarize the approximation error per model
summary_list <- lapply(split(lambda_check, lambda_check$model), function(d) {
  error <- d$lambda.exact - d$lambda.approx
  viable <- d$lambda.exact > 0
  data.frame(
    model = d$model[1],
    n = nrow(d),
    n.viable = sum(viable),
    max.abs.error = max(abs(error)),
    rmse = sqrt(mean(error^2)),
    median.rel.error.viable = median(abs(error[viable]) / d$lambda.exact[viable]),
    r.squared = cor(d$lambda.approx, d$lambda.exact)^2
  )
})
approx_summary <- do.call(rbind, summary_list)
print(approx_summary)

write.table(approx_summary, "./data/lambda_approximation_summary.txt", sep = "\t", row.names = F, quote = F)

##########

pdf("./output/Fig_S10_lambda_approximation_validation.pdf", width = 2.5, height = 2.5)
#### Plot approximated versus exact lambda ----

site_colors <- c(california = "#1f78b4", brazil = "#33a02c", yemen = "#ff7f00")

for (model_type in c("composite", "separate")) {

  model_dat <- lambda_check[lambda_check$model == model_type, ]

  # R2 treating x as the exact prediction of y (slope 1, intercept 0), all points
  r2_pred <- 1 - sum((model_dat$lambda.exact^24 - model_dat$lambda.approx^24)^2) /
    sum((model_dat$lambda.exact^24 - mean(model_dat$lambda.exact^24))^2)

  # Thin the composite point cloud for plotting
  if (nrow(model_dat) > 20000) model_dat <- model_dat[sample(1:nrow(model_dat), 20000), ]

  # Daily scale (^24) accentuates differences near lambda = 1, where they matter most
  plot(model_dat$lambda.approx^24, model_dat$lambda.exact^24,
       col = paste0(site_colors[model_dat$location], "30"),
       pch = 16, cex = 0.5,
       xlim = c(0, 1.25), ylim = c(0, 1.25),
       xlab = "Approximated lambda (daily)",
       ylab = "Euler-Lotka lambda (daily)",
       main = paste0(model_type, " model"))
  abline(0, 1, lty = 2)
  text(0, 1.25, bquote(italic(R)^2 == .(format(round(r2_pred, 7), nsmall = 7))), adj = c(0, 1))
}

##########
dev.off()

#### Sensitivity 2: timing of infertility risk accumulation (Eq. 9) ----

# Reallocate all infertility risk to the first or last 10% of development
# (rescaled so constant-temperature outcomes are unchanged) and compare with
# the uniform accumulation assumed by Eq. 9, for the same 50 draws

# Mirrors separate_fitness_function() (functions.R), adding window-concentrated fertility
separate_lambda_fert_window <- function(temperatures, origin) {

  n_temps <- length(temperatures)
  main_cohorts <- 1:365
  main_start_hours <- (364 + main_cohorts) * 24 + 1
  if (max(main_start_hours) > n_temps) stop("Temperature vector too short.")

  all_dev_rates <- development_rate(temperatures, origin)
  all_growth_rates <- growth_rate(temperatures, origin)
  all_surv_rates <- survival_rate(temperatures, origin)
  all_fert_rates <- fertility_rate(temperatures, origin)

  C_lethal_count <- c(0, cumsum(all_surv_rates <= 0))
  C_infertile_count <- c(0, cumsum(all_fert_rates <= 0))

  all_surv_rates[all_surv_rates <= 0] <- 1e-16
  all_fert_rates[all_fert_rates <= 0] <- 1e-16

  all_surv_rates_log <- log(all_surv_rates)
  all_fert_rates_log <- log(all_fert_rates)
  ovi_temp_factors <- fecundity_to_oviposition_rate(1, temperatures)

  C_dev <- c(0, cumsum(all_dev_rates))
  C_growth <- c(0, cumsum(all_growth_rates))
  C_surv_log <- c(0, cumsum(all_surv_rates_log))
  C_fert_log <- c(0, cumsum(all_fert_rates_log))

  dev_targets <- 1 + C_dev[main_start_hours]
  natural_end_hours <- findInterval(dev_targets, C_dev)
  raw_duration <- natural_end_hours - main_start_hours + 1
  max_duration_hours <- 250 * 24
  simulated_dev_time_hours <- ifelse(is.na(raw_duration) | raw_duration > max_duration_hours,
                                     max_duration_hours,
                                     raw_duration)
  dev_end_hours <- main_start_hours + simulated_dev_time_hours - 1
  dev_end_hours[dev_end_hours >= n_temps] <- n_temps

  simulated_mass <- C_growth[dev_end_hours + 1] - C_growth[main_start_hours]
  simulated_mass[simulated_mass < 0.0005] <- 0.0005
  simulated_mass[simulated_mass > 0.0035] <- 0.0035

  total_survival <- exp(C_surv_log[dev_end_hours + 1] - C_surv_log[main_start_hours])
  lethal_hits <- C_lethal_count[dev_end_hours + 1] - C_lethal_count[main_start_hours]
  total_survival[lethal_hits > 0] <- 0

  # Uniform risk accumulation (Eq. 9, as in separate_fitness_function)
  total_fertility <- exp(C_fert_log[dev_end_hours + 1] - C_fert_log[main_start_hours])
  total_fertility[is.na(total_fertility)] <- 0
  infertile_hits <- C_infertile_count[dev_end_hours + 1] - C_infertile_count[main_start_hours]
  total_fertility[infertile_hits > 0] <- 0

  # Hours at which cumulative development crosses 10% and 90%
  dev_total <- C_dev[dev_end_hours + 1] - C_dev[main_start_hours]
  h10 <- pmin(findInterval(C_dev[main_start_hours] + 0.1, C_dev), dev_end_hours)
  h90 <- pmin(findInterval(C_dev[main_start_hours] + 0.9, C_dev), dev_end_hours)
  reached_late <- dev_total >= 0.9

  # Risk in first 10%: log-hazard rescaled by dev_total / dev_window
  dev_early <- C_dev[h10 + 1] - C_dev[main_start_hours]
  log_f_early <- ifelse(dev_early > 0,
                        (dev_total / dev_early) * (C_fert_log[h10 + 1] - C_fert_log[main_start_hours]),
                        0)
  fert_early <- exp(log_f_early)
  fert_early[is.na(fert_early)] <- 0
  hits_early <- C_infertile_count[h10 + 1] - C_infertile_count[main_start_hours]
  fert_early[hits_early > 0] <- 0

  # Risk in last 10%: no risk accrues if development never reaches 90%
  dev_late <- C_dev[dev_end_hours + 1] - C_dev[h90]
  log_f_late <- ifelse(reached_late & dev_late > 0,
                       (dev_total / dev_late) * (C_fert_log[dev_end_hours + 1] - C_fert_log[h90]),
                       0)
  fert_late <- exp(log_f_late)
  fert_late[is.na(fert_late)] <- 0
  hits_late <- C_infertile_count[dev_end_hours + 1] - C_infertile_count[h90]
  fert_late[reached_late & hits_late > 0] <- 0

  lambda_uniform <- numeric(365)
  lambda_early <- numeric(365)
  lambda_late <- numeric(365)

  for (i in 1:365) {

    total_egg_pool <- mass_to_fecundity(simulated_mass[i], origin) / 2

    adult_start_hour <- dev_end_hours[i] + 1
    eggs_laid_so_far <- 0
    hour_counter <- 0
    hourly_offspring_vector <- c()

    while (eggs_laid_so_far < total_egg_pool) {
      current_hour <- adult_start_hour + hour_counter
      if (current_hour > n_temps) break

      eggs_this_hour <- total_egg_pool * ovi_temp_factors[current_hour]
      # Clamp eggs to remaining pool
      eggs_remaining <- total_egg_pool - eggs_laid_so_far
      if(eggs_this_hour > eggs_remaining) eggs_this_hour <- eggs_remaining

      hourly_offspring_vector <- c(hourly_offspring_vector, eggs_this_hour)
      eggs_laid_so_far <- eggs_laid_so_far + eggs_this_hour
      hour_counter <- hour_counter + 1

      if(hour_counter > 200 * 24) break
    }

    age_of_mother_hours <- simulated_dev_time_hours[i] + (1:length(hourly_offspring_vector))
    mx <- hourly_offspring_vector
    x_hours <- age_of_mother_hours
    lx <- total_survival[i] * total_fertility[i]

    R0 <- sum(mx) * lx
    T_gen_hours <- sum(x_hours * mx) / sum(mx)
    lambda_uniform[i] <- R0 ^ (1 / (T_gen_hours))

    lambda_early[i] <- (sum(mx) * total_survival[i] * fert_early[i]) ^ (1 / (T_gen_hours))
    lambda_late[i] <- (sum(mx) * total_survival[i] * fert_late[i]) ^ (1 / (T_gen_hours))
  }

  lambda_uniform[is.na(lambda_uniform)] <- 0
  lambda_early[is.na(lambda_early)] <- 0
  lambda_late[is.na(lambda_late)] <- 0

  return(list(lambda.uniform = lambda_uniform,
              lambda.early10 = lambda_early,
              lambda.late10 = lambda_late,
              fert.uniform = total_fertility,
              fert.early10 = fert_early,
              fert.late10 = fert_late))
}

results_list2 <- list()

for (draw in 1:50) {

  location_data <- subset(temps, unique.ID == sampled_sites[draw])
  location_name <- location_data$location[1]
  warmed_temperatures <- location_data$temperature + sampled_warming[draw]

  res <- separate_lambda_fert_window(warmed_temperatures, location_name)

  # Cross-check: the uniform variant must reproduce the pipeline output
  stopifnot(isTRUE(all.equal(res$lambda.uniform,
                             separate_fitness_function(warmed_temperatures, location_name)$lambda.hourly)))

  results_list2[[draw]] <- data.frame(
    unique.ID = sampled_sites[draw],
    location = location_name,
    warming = sampled_warming[draw],
    lambda.uniform = res$lambda.uniform,
    lambda.early10 = res$lambda.early10,
    lambda.late10 = res$lambda.late10,
    fert.uniform = res$fert.uniform,
    fert.early10 = res$fert.early10,
    fert.late10 = res$fert.late10
  )

  cat(paste("  - Completed draw", draw, "of 50:", sampled_sites[draw], "+", sampled_warming[draw], "C\n"))
}

fert_window_check <- do.call(rbind, results_list2)

# Summarize the deviation from the uniform assumption per variant
summary_list2 <- lapply(c("early10", "late10"), function(w) {
  lambda_diff <- fert_window_check[[paste0("lambda.", w)]] - fert_window_check$lambda.uniform
  fert_diff_pct <- (fert_window_check[[paste0("fert.", w)]] - fert_window_check$fert.uniform) * 100
  data.frame(
    variant = w,
    n = length(lambda_diff),
    max.abs.lambda.diff = max(abs(lambda_diff)),
    rmse.lambda = sqrt(mean(lambda_diff^2)),
    r.squared.lambda = cor(fert_window_check$lambda.uniform, fert_window_check[[paste0("lambda.", w)]])^2,
    median.abs.fert.diff.pct = median(abs(fert_diff_pct)),
    max.abs.fert.diff.pct = max(abs(fert_diff_pct)),
    n.sterile.uniform.only = sum(fert_window_check$lambda.uniform == 0 & fert_window_check[[paste0("lambda.", w)]] > 0),
    n.sterile.variant.only = sum(fert_window_check$lambda.uniform > 0 & fert_window_check[[paste0("lambda.", w)]] == 0)
  )
})
fert_window_summary <- do.call(rbind, summary_list2)
print(fert_window_summary)

write.table(fert_window_summary, "./data/fertility_window_summary.txt", sep = "\t", row.names = F, quote = F)

##########

pdf("./output/Fig_S10_fertility_window_sensitivity.pdf", width = 2.5, height = 2.5)
#### Plot uniform versus window-concentrated infertility risk ----

# Daily-scale lambda: uniform (Eq. 9) on x, window-concentrated variants on y
for (variant in c("lambda.early10", "lambda.late10")) {

  window_label <- ifelse(variant == "lambda.early10", "first", "last")

  # R2 treating x as the exact prediction of y (slope 1, intercept 0)
  r2_pred <- 1 - sum((fert_window_check[[variant]]^24 - fert_window_check$lambda.uniform^24)^2) /
    sum((fert_window_check[[variant]]^24 - mean(fert_window_check[[variant]]^24))^2)

  plot(fert_window_check$lambda.uniform^24, fert_window_check[[variant]]^24,
       col = paste0(site_colors[fert_window_check$location], "30"),
       pch = 16, cex = 0.5,
       xlim = c(0, 1.25), ylim = c(0, 1.25),
       xlab = "Lambda, uniform risk (daily)",
       ylab = paste0("Lambda, risk in ", window_label, " 10% (daily)"),
       main = paste0("Risk in ", window_label, " 10%"))
  abline(0, 1, lty = 2)
  text(0, 1.25, bquote(italic(R)^2 == .(format(round(r2_pred, 4), nsmall = 4))), adj = c(0, 1))
}

##########
dev.off()

#### Sensitivity 3: temperatures within the Tmax uncertainty range ----

# Share of hours falling within the 90% credible interval (HPDI, cf.
# parameter_estimates_*.txt) for Tmax of the development and growth rate
# TPCs, where it is uncertain whether temperatures exceed Tmax, under
# present-day conditions and +5 C warming

site_temps <- split(temps$temperature, temps$unique.ID)[eligible_sites]
site_locations <- temps$location[match(eligible_sites, temps$unique.ID)]

origin_codes <- c(california = "ca", brazil = "bra", yemen = "yem")
tmax_intervals <- list(development = dev_parameters[dev_parameters$parameter == "Tmax", ],
                       growth = growth_parameters[growth_parameters$parameter == "Tmax", ])

# Per-site share of hours within the origin-specific Tmax interval
share_within_interval <- function(trait, warming) {
  ints <- tmax_intervals[[trait]]
  vapply(seq_along(site_temps), function(s) {
    lo <- ints$hpdi.lower[ints$origin == origin_codes[site_locations[s]]]
    hi <- ints$hpdi.upper[ints$origin == origin_codes[site_locations[s]]]
    warmed <- site_temps[[s]] + warming
    mean(warmed >= lo & warmed <= hi) * 100
  }, numeric(1))
}

bar_specs <- expand.grid(trait = c("development", "growth"), warming = c(0, 5))
site_shares <- lapply(seq_len(nrow(bar_specs)), function(i) {
  share_within_interval(as.character(bar_specs$trait[i]), bar_specs$warming[i])
})

exposure_summary <- data.frame(bar_specs,
                               mean.pct.hours = sapply(site_shares, mean),
                               median.pct.hours = sapply(site_shares, median),
                               max.pct.hours = sapply(site_shares, max))
print(exposure_summary)

write.table(exposure_summary, "./data/extrapolation_exposure_summary.txt", sep = "\t", row.names = F, quote = F)

##########

pdf("./output/Fig_S11_extrapolation_uncertainty.pdf", width = 2.5, height = 2.5)
#### Plot the share of hours within the Tmax uncertainty range ----

# Bars are means across sites; points are individual sites, colored by region
bar_pos <- barplot(exposure_summary$mean.pct.hours,
                   space = c(0.2, 0.2, 0.8, 0.2),
                   col = c("#666666", "#bbbbbb"),
                   ylim = c(0, max(unlist(site_shares))),
                   ylab = "Hours within Tmax 90% CrI (%)")
mtext(rep(c("dev.", "growth"), 2), side = 1, line = 0.3, at = bar_pos, cex = 0.55)
mtext(c("present-day", "+5°C"), side = 1, line = 1.5,
      at = c(mean(bar_pos[1:2]), mean(bar_pos[3:4])), cex = 0.8)

for (i in seq_along(site_shares)) {
  points(jitter(rep(bar_pos[i], length(site_shares[[i]])), amount = 0.2), site_shares[[i]],
         pch = 16, cex = 0.4, col = paste0(site_colors[site_locations], "30"))
}

##########
dev.off()
