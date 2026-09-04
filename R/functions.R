# Read thermal response parameters
LRS_parameters <- read.delim("./data/parameter_estimates_LRS.txt")
growth_parameters <- read.delim("./data/parameter_estimates_growth_rate.txt")
dev_parameters <- read.delim("./data/parameter_estimates_dev_rate.txt")

# Create a consolidated list of fecundity parameters
fec_params <- list(
  yemen = list(
    Topt = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "ToptFec"],
    fecmax = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "fecmax"],
    Tb = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "Tb"],
    K = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "K"]),
  california = list(
    Topt = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "ToptFec"],
    fecmax = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "fecmax"],
    Tb = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "Tb"],
    K = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "K"]),
  brazil = list(
    Topt = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "ToptFec"],
    fecmax = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "fecmax"],
    Tb = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "Tb"],
    K = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "K"])
)

# Create a consolidated list of growth parameters
growth_params <- list(
  yemen = list(
    Tmin = growth_parameters$MAP[growth_parameters$origin == "yem" & growth_parameters$parameter == "Tmin"],
    Topt = growth_parameters$MAP[growth_parameters$origin == "yem" & growth_parameters$parameter == "Topt"],
    Tmax = growth_parameters$MAP[growth_parameters$origin == "yem" & growth_parameters$parameter == "Tmax"],
    rmax = growth_parameters$MAP[growth_parameters$origin == "yem" & growth_parameters$parameter == "rmax"]
  ),
  california = list(
    Tmin = growth_parameters$MAP[growth_parameters$origin == "ca" & growth_parameters$parameter == "Tmin"],
    Topt = growth_parameters$MAP[growth_parameters$origin == "ca" & growth_parameters$parameter == "Topt"],
    Tmax = growth_parameters$MAP[growth_parameters$origin == "ca" & growth_parameters$parameter == "Tmax"],
    rmax = growth_parameters$MAP[growth_parameters$origin == "ca" & growth_parameters$parameter == "rmax"]
  ),
  brazil = list(
    Tmin = growth_parameters$MAP[growth_parameters$origin == "bra" & growth_parameters$parameter == "Tmin"],
    Topt = growth_parameters$MAP[growth_parameters$origin == "bra" & growth_parameters$parameter == "Topt"],
    Tmax = growth_parameters$MAP[growth_parameters$origin == "bra" & growth_parameters$parameter == "Tmax"],
    rmax = growth_parameters$MAP[growth_parameters$origin == "bra" & growth_parameters$parameter == "rmax"]
  )
)

# Create a consolidated list of development parameters
dev_params <- list(
  yemen = list(
    Tmin = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "Tmin"],
    Topt = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "Topt"],
    Tmax = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "Tmax"],
    rmax = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "rmax"]
  ),
  california = list(
    Tmin = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "Tmin"],
    Topt = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "Topt"],
    Tmax = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "Tmax"],
    rmax = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "rmax"]
  ),
  brazil = list(
    Tmin = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "Tmin"],
    Topt = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "Topt"],
    Tmax = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "Tmax"],
    rmax = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "rmax"]
  )
)

# Create a consolidated list of fertility parameters
fert_params <- list(
  yemen = list(
    Tmin = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "Tmin"],
    Topt = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "Topt"],
    Tmax = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "Tmax"],
    rmax = dev_parameters$MAP[dev_parameters$origin == "yem" & dev_parameters$parameter == "rmax"],
    a    = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "a"],
    h    = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "ToptViab"],
    k    = LRS_parameters$MAP[LRS_parameters$origin == "yem" & LRS_parameters$parameter == "k"]
  ),
  california = list(
    Tmin = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "Tmin"],
    Topt = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "Topt"],
    Tmax = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "Tmax"],
    rmax = dev_parameters$MAP[dev_parameters$origin == "ca" & dev_parameters$parameter == "rmax"],
    a    = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "a"],
    h    = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "ToptViab"],
    k    = LRS_parameters$MAP[LRS_parameters$origin == "ca" & LRS_parameters$parameter == "k"]
  ),
  brazil = list(
    Tmin = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "Tmin"],
    Topt = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "Topt"],
    Tmax = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "Tmax"],
    rmax = dev_parameters$MAP[dev_parameters$origin == "bra" & dev_parameters$parameter == "rmax"],
    a    = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "a"],
    h    = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "ToptViab"],
    k    = LRS_parameters$MAP[LRS_parameters$origin == "bra" & LRS_parameters$parameter == "k"]
  )
)


# Geometric mean
geo_mean <- function(x) {
  exp(mean(log(x)))
}

# The inverse logit function 
invlogit <- function(x) (1 / (1 + exp(-x)))

# LRF TPC function for curve fitting
LRF <- function(temp, Tmin, Topt, Tmax, rmax) {
  # Create a logical index for temperatures within the valid range
  valid_idx <- temp >= Tmin & temp <= Tmax
  
  # Initialize the result vector with zeros
  res <- numeric(length(temp))
  
  # Only perform calculations for valid temperatures to avoid errors and save time
  if (!any(valid_idx)) return(res)
  
  # Subset the temperature vector for calculation
  temp_valid <- temp[valid_idx]
  
  # Convert temperatures to Kelvin
  temp_k <- temp_valid + 273.15
  Tmin_k <- Tmin + 273.15
  Tmax_k <- Tmax + 273.15
  Topt_k <- Topt + 273.15
  
  # Perform the calculation on the subset of valid temperatures
  numerator <- rmax * (temp_k - Tmax_k) * (temp_k - Tmin_k)^2
  
  term1 <- (Topt_k - Tmin_k)
  term2 <- term1 * (temp_k - Topt_k)
  term3 <- (Topt_k - Tmax_k) * (Topt_k + Tmin_k - 2 * temp_k)
  denominator <- term1 * (term2 - term3)
  
  # Use a small tolerance for floating point comparison
  calculated_vals <- ifelse(abs(denominator) < 1e-9, 0, numerator / denominator)
  
  # Place the calculated values back into the result vector, ensuring no negative rates
  res[valid_idx] <- pmax(0, calculated_vals)
  
  return(res)
}

# Non-viability function for curve fitting
nonviability <- function(temp, a, Topt, k) {
  a * (temp - Topt)^2 + k
}

# Fecundity function for curve fitting
fecundity <- function(temp, Topt, fecmax, Tb, K) {
  20 ^ (-abs((temp - Topt) / (Tb / 2))^(K + 1)) * fecmax
}

# Specify base fecundity function
base_fecundity <- function(temperature, origin) {
  params <- fec_params[[origin]]
  fec <- fecundity(temp = temperature,
                   Topt = params$Topt,
                   fecmax = params$fecmax,
                   Tb = params$Tb,
                   K = params$K)
  return(fec)
}

# Fertility rate functions
# General 
fert_rate <- function(temp, Tmin, Topt, Tmax, rmax, a, h, k) {
  # Clamp temperatures for the development rate calculation 
  dev_temp <- pmin(pmax(temp, 17), 37)
  
  # Use the LRF function
  dev_rate <- LRF(dev_temp, Tmin, Topt, Tmax, rmax)
  
  prop_viable <- 1 - invlogit(nonviability(temp, a, h, k))
  
  return(prop_viable^dev_rate)
}

# Fixed
fertility_rate <- function(temperature, origin) {
  params <- fert_params[[origin]]

  rate <- fert_rate(temp = temperature, 
                    Tmin = params$Tmin, Topt = params$Topt, Tmax = params$Tmax, rmax = params$rmax,
                    a = params$a, h = params$h, k = params$k)
  return(rate)
}

# Development rate rate function (fixed)
development_rate <- function(temperature, origin) {
  
  params <- dev_params[[origin]]
  
  rate <- LRF(temp = temperature, 
              Tmin = params$Tmin, 
              Topt = params$Topt, 
              Tmax = params$Tmax, 
              rmax = params$rmax)
  
  return(rate)
}

# Growth rate rate function (fixed)
growth_rate <- function(temperature, origin) {
  
  params <- growth_params[[origin]]
  
  rate <- LRF(temp = temperature, 
              Tmin = params$Tmin, 
              Topt = params$Topt, 
              Tmax = params$Tmax, 
              rmax = params$rmax)
  
  return(rate)
}


# Survival rate function (fixed)
survival_rate <- function(temperature, origin){
  
  # Calculate devrate and mass using pmin/pmax for vector safety
  bounded_temp <- pmin(pmax(temperature, 17), 37)
  devrate <- development_rate(bounded_temp, origin)
  
  # Calculate mass
  raw_mass <- growth_rate(bounded_temp, origin) / devrate
  mass <- pmax(raw_mass, 0.0005)
  
  max_fecundity <- mass_to_fecundity(mass, origin)
  surviving_proportion <- base_fecundity(temperature, origin) / max_fecundity
  surviving_proportion <- pmin(surviving_proportion, 1)
  
  survrate <- surviving_proportion^devrate
  
  # Handle NaNs (Vectorized replacement)
  survrate[is.nan(survrate)] <- 0
  
  return(survrate)
}


# Mass to max fecundity function
# conversion_fec <- lm(log(adult.offspring) ~ 0 + log(average.mass) + origin, data = subset(dat, temperature == 29))
# coef(conversion_fec)
mass_to_fecundity <- function(mass, origin){
  if(origin == "california") { intercept = 6.9582461 }
  if(origin == "brazil") {intercept = 7.1822847 }
  if(origin == "yemen") { intercept = 7.2572643 }
  return(exp(intercept + 0.4440254 * log(mass)))
}

# Fecundity to oviposition rate function
# conversion_ovi <- lm(log(eggs.1h / LRS) ~ temperature, data = ovi_conversion)
# coef(conversion_ovi)
fecundity_to_oviposition_rate <- function(fecundity, temperature){
  return(fecundity * exp(-5.49991935 +  0.06742549 * temperature))
}

# Create the composite function
composite_fitness_function <- function(temperature, origin, IGR = T, return_components = FALSE) {
  
  # 1. Get hourly physiological rates
  dev_r <- development_rate(temperature = temperature, origin = origin)
  growth_r <- growth_rate(temperature = temperature, origin = origin)
  surv_r <- survival_rate(temperature = temperature, origin = origin)
  fert_r <- fertility_rate(temperature = temperature, origin = origin)
  
  # 2. Calculate development time (hours)
  raw_dev_time <- 1 / dev_r
  T_dev_hours <- ceiling(raw_dev_time)
  T_dev_hours <- ifelse(!is.na(T_dev_hours) & T_dev_hours > 250 * 24, 250 * 24, T_dev_hours)  
  
  # 3. Calculate final mass (g)
  mass_ind <- growth_r * T_dev_hours
  mass_ind <- ifelse(!is.na(mass_ind) & mass_ind < 0.0005, 0.0005, mass_ind) 
  mass_ind <- ifelse(!is.na(mass_ind) & mass_ind > 0.0035, 0.0035, mass_ind)
  
  # 4. Calculate survival and fertility probabilities (Total for dev period)
  survival_total <- surv_r ^ T_dev_hours
  survival_total[is.na(survival_total)] <- 0 
  fertility_total <- fert_r ^ T_dev_hours
  fertility_total[is.na(fertility_total)] <- 0 
  
  # 5. Reproductive traits
  max_fec <- mass_to_fecundity(mass = mass_ind, origin = origin) / 2 
  ovi_r <- fecundity_to_oviposition_rate(fecundity = max_fec, temperature = temperature)
  
  raw_ovi_time <- max_fec / ovi_r
  T_adult_hours <- ceiling(raw_ovi_time)
  
  # 6. Calculate Fitness (Hourly Lambda)
  R0 <- survival_total * fertility_total * max_fec
  
  # Oviposition moments
  ovi_fraction <- fecundity_to_oviposition_rate(1, temperature)
  raw_duration <- 1 / ovi_fraction
  k_full_hours <- floor(raw_duration)
  
  remainder_weight <- 1 - (k_full_hours * ovi_fraction)
  remainder_weight <- ifelse(remainder_weight < 0, 0, remainder_weight)
  
  moment_full <- ovi_fraction * (k_full_hours * (k_full_hours + 1) / 2)
  moment_remainder <- remainder_weight * (k_full_hours + 1)
  mean_ovi_time <- moment_full + moment_remainder
  T_gen_hours <- T_dev_hours + mean_ovi_time
  
  # Hourly Lambda Calculation
  fitness_lambda_hourly <- R0^(1 / T_gen_hours)
  fitness_lambda_hourly[is.na(fitness_lambda_hourly)] <- 0
  
  
  final_fitness <- if(IGR) fitness_lambda_hourly else c(max_fec * survival_total * fertility_total * 2/ T_dev_hours)
  
  if (return_components) {
    
    return(list(
      lambda.hourly = final_fitness,
      dev.time = T_dev_hours / 24,
      mass = mass_ind,
      fertility.daily = fert_r^24,
      R0 = R0,
      ovi.time.days = T_adult_hours / 24,
      survival.daily = surv_r^24
      
    ))
  } else {
    return(final_fitness)
  }
}

# Create the separate components function
separate_fitness_function <- function(temperatures, origin) {
  
  n_temps <- length(temperatures)
  main_cohorts <- 1:365
  main_start_hours <- (364 + main_cohorts) * 24 + 1 
  if (max(main_start_hours) > n_temps) stop("Temperature vector too short.")
  
  # 1. Pre-calculate all hourly rates 
  all_dev_rates <- development_rate(temperatures, origin)
  all_growth_rates <- growth_rate(temperatures, origin)
  all_surv_rates <- survival_rate(temperatures, origin)
  all_fert_rates <- fertility_rate(temperatures, origin)
  
  # Create a "death counter"
  C_lethal_count <- c(0, cumsum(all_surv_rates <= 0))
  C_infertile_count <- c(0, cumsum(all_fert_rates <= 0))
  
  # Replace 0s with small number
  all_surv_rates[all_surv_rates <= 0] <- 1e-16
  all_fert_rates[all_fert_rates <= 0] <- 1e-16
  
  all_surv_rates_log <- log(all_surv_rates)
  all_fert_rates_log <- log(all_fert_rates)
  ovi_temp_factors <- fecundity_to_oviposition_rate(1, temperatures)
  
  C_dev <- c(0, cumsum(all_dev_rates))
  C_growth <- c(0, cumsum(all_growth_rates))
  C_surv_log <- c(0, cumsum(all_surv_rates_log))
  C_fert_log <- c(0, cumsum(all_fert_rates_log))
  
  # 2. Determine development end times
  dev_targets <- 1 + C_dev[main_start_hours]
  natural_end_hours <- findInterval(dev_targets, C_dev)
  raw_duration <- natural_end_hours - main_start_hours + 1
  max_duration_hours <- 250 * 24
  simulated_dev_time_hours <- ifelse(is.na(raw_duration) | raw_duration > max_duration_hours,
                                     max_duration_hours,
                                     raw_duration)
  dev_end_hours <- main_start_hours + simulated_dev_time_hours - 1
  dev_end_hours[dev_end_hours >= n_temps] <- n_temps
  
  # 3. Calculate accumulated traits
  simulated_mass <- C_growth[dev_end_hours + 1] - C_growth[main_start_hours]
  simulated_mass[simulated_mass < 0.0005] <- 0.0005
  simulated_mass[simulated_mass > 0.0035] <- 0.0035
  
  # Total Survival/Fertility over development
  total_survival <- exp(C_surv_log[dev_end_hours + 1] - C_surv_log[main_start_hours])
  lethal_hits <- C_lethal_count[dev_end_hours + 1] - C_lethal_count[main_start_hours]
  total_survival[lethal_hits > 0] <- 0
  
  total_fertility <- exp(C_fert_log[dev_end_hours + 1] - C_fert_log[main_start_hours])
  total_fertility[is.na(total_fertility)] <- 0
  infertile_hits <- C_infertile_count[dev_end_hours + 1] - C_infertile_count[main_start_hours]
  total_fertility[infertile_hits > 0] <- 0
  
  
  # Calculate Daily Probabilities: (Total)^(24/T_hours)
  # This normalizes risk to a 24hr period
  daily_survival_prob <- total_survival ^ (24 / simulated_dev_time_hours)
  daily_fertility_prob <- total_fertility ^ (24 / simulated_dev_time_hours)
  
  simulated_lambda <- numeric(365)
  simulated_offspring <- numeric(365)
  simulated_ovi_time_hours <- rep(NA, 365)
  simulated_max_fecundity <- rep(NA, 365)
  
  # 4. Loop through daily cohorts
  for (i in 1:365) {
    
    total_egg_pool <- mass_to_fecundity(simulated_mass[i], origin) / 2
    simulated_max_fecundity[i] <- total_egg_pool
    
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
    
    simulated_ovi_time_hours[i] <- hour_counter
    
    age_of_mother_hours <- simulated_dev_time_hours[i] + (1:length(hourly_offspring_vector))
    mx <- hourly_offspring_vector
    x_hours <- age_of_mother_hours
    lx <- total_survival[i] * total_fertility[i]
    simulated_offspring[i] <- sum(hourly_offspring_vector, na.rm = TRUE) * lx 
    
    R0 <- sum(mx) * lx
    T_gen_hours <- sum(x_hours * mx) / sum(mx)
    simulated_lambda[i] <- R0 ^ (1 / (T_gen_hours))
    
    
  }
  
  simulated_lambda[is.na(simulated_lambda)] <- 0
  
  return(list(
    lambda.hourly =  simulated_lambda,
    dev.time = simulated_dev_time_hours / 24,
    mass = simulated_mass, 
    fertility.daily = daily_fertility_prob, # Returns daily prob
    R0 = simulated_offspring,
    ovi.time.days = simulated_ovi_time_hours / 24,
    survival.daily = daily_survival_prob   # Returns daily prob
  ))
}

# Main wrapper (Summary and Plotting)
calculate_fitness <- function(temperatures, origin, plot = FALSE, return_hourly = FALSE, return_all = FALSE) {
  
  separate_results_list <- list()
  composite_results_list <- list() 
  get_components <- plot | return_all
  
  year_hours <- 365 * 24
  temps_year2 <- if(length(temperatures) >= (2 * year_hours)) temperatures[(year_hours + 1):(2 * year_hours)] else temperatures 
  
  for (ori_name in origin) {
    # Run separate
    separate_results_list[[ori_name]] <- separate_fitness_function(temperatures, ori_name)
    # Run composite
    composite_results_list[[ori_name]] <- composite_fitness_function(temps_year2, ori_name, return_components = get_components)
  }
  
  # Summary DataFrame
  summary_list <- lapply(names(separate_results_list), function(name) {
    sep_res <- separate_results_list[[name]]
    comp_res <- composite_results_list[[name]] 
    
    
    if (get_components) {
      # Use the pre-calculated daily metrics from composite
      data.frame(
        origin = name,
        lambda.separate.arith = mean(sep_res$lambda.hourly, na.rm = TRUE),
        lambda.separate.geo = geo_mean(sep_res$lambda.hourly),
        lambda.composite.arith = mean(comp_res$lambda.hourly, na.rm = TRUE), 
        lambda.composite.geo = geo_mean(comp_res$lambda.hourly),          
        mass.separate = mean(sep_res$mass, na.rm = TRUE),
        mass.composite = mean(comp_res$mass, na.rm = TRUE),
        dev.time.separate = mean(sep_res$dev.time, na.rm = TRUE),
        dev.time.composite = mean(comp_res$dev.time, na.rm = TRUE),
        R0.separate = mean(sep_res$R0, na.rm = TRUE),
        R0.composite = mean(comp_res$R0, na.rm = TRUE),
        ovi.time.separate = mean(sep_res$ovi.time.days, na.rm = TRUE),
        ovi.time.composite = mean(comp_res$ovi.time.days, na.rm = TRUE),
        fertility.separate = mean(sep_res$fertility.daily, na.rm = TRUE),
        fertility.composite = mean(comp_res$fertility.daily, na.rm = TRUE),
        survival.separate = mean(sep_res$survival.daily, na.rm = TRUE),
        survival.composite = mean(comp_res$survival.daily, na.rm = TRUE)
      )
    } else {
      
      data.frame(
        origin = name,
        lambda.separate.ari = mean(sep_res$lambda.hourly, na.rm = TRUE),
        lambda.separate.geo = geo_mean(sep_res$lambda.hourly),
        lambda.composite.ari = mean(comp_res, na.rm = TRUE),
        lambda.composite.geo = geo_mean(comp_res)
      )
    }
  })
  results_df <- do.call(rbind, summary_list)
  
  # Plotting logic
  if (plot) {
    par(mfrow = c(3, 3), mar = c(4.1, 4.1, 3.1, 1.1), oma = c(0, 0, 2, 0))
    
    # Define plot targets. 
    plot_vars <- list(
      list(id = "dev.time", title = "Development time (days)"),
      list(id = "mass", title = "Final mass (g)"),
      list(id = "fertility.daily", title = "Daily Fertility Prob."),
      list(id = "R0", title = "R0"),
      list(id = "survival.daily", title = "Daily Survival Prob."),
      list(id = "ovi.time.days", title = "Oviposition time (days)"),
      list(id = "lambda.hourly", title = "Daily Lambda (Arith Mean)", sep_id = "lambda.hourly"),
      list(id = "lambda.hourly", title = "Daily Lambda (Geo Mean)", sep_id = "lambda.hourly")
    )
    
    for (p_var in plot_vars) {
      sep_key <- if(!is.null(p_var$sep_id)) p_var$sep_id else p_var$id
      sep_data <- lapply(separate_results_list, `[[`, sep_key)
      comp_data <- lapply(composite_results_list, `[[`, p_var$id)
      
      valid_vals <- c(unlist(sep_data), unlist(comp_data))
      valid_vals <- valid_vals[is.finite(valid_vals)]
      ylim <- if(length(valid_vals) == 0) c(0,1) else c(min(valid_vals, na.rm=TRUE) * 0.9, max(valid_vals, na.rm=TRUE) * 1.1)
      
      plot(NA, xlim = c(1, 365), ylim = ylim, xlab = "Day of Year", ylab = "", main = p_var$title)
      
      for (name in names(composite_results_list)) {
        c_dat <- comp_data[[name]]
        lines((1:(365*24))/24, c_dat, col = "#FF00FF20", lwd = 1)
      }
      for (name in names(separate_results_list)) {
        lines(1:365, sep_data[[name]], col = "black", lwd = 2)
      }
      
    }
    
    # Temp plot
    daily_matrix <- matrix(temps_year2, nrow = 24)
    daily_means <- colMeans(daily_matrix, na.rm = TRUE)
    daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
    daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)
    
    plot(NA, 
         xlim = c(1, 365), 
         ylim = range(c(daily_mins, daily_maxs), na.rm = TRUE), 
         xlab = "Day of Year", 
         ylab = "Temp (°C)", 
         main = "Temperature Regime (Year 2)")
    
    if(origin[1] == "california") {color <- "#1f78b4"} else if(origin[1] == "brazil") {color <- "#33a02c"} else {color <- "#000000"}
    
    polygon(x = c(1:365, 365:1),
            y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
            col = paste0(color, "40"),
            border = NA)
    lines(1:365, daily_means[1:365], col = color, lwd = 2)      
    
    mtext("Separate components (black) vs composite function red)", outer = TRUE, cex = 1.5)
    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  }
  
  if (return_hourly) {
    daily_results_list <- lapply(names(separate_results_list), function(name) {
      comp_res <- composite_results_list[[name]]
      data.frame(
        origin = name, hour = 1:(365*24),
        lambda.separate = rep(separate_results_list[[name]]$lambda.hourly, each = 24),
        lambda.composite = comp_res[1:(365*24)]
      )
    })
    return(list(summary = results_df, hourly = do.call(rbind, daily_results_list)))
  } else {
    return(results_df)
  }
}


if (paste0("package:brms") %in% search()) {
stan_funs <- "
  /**
   * LRF Thermal Performance Curve (TPC) function.
   */
  real LRF(real temp, real Tmin, real Topt, real Tmax, real rmax) {
    // Return 0 if temperature is outside the valid range [Tmin, Tmax]
    if (temp < Tmin || temp > Tmax) {
      return 0.0;
    }

    // Convert temperatures to Kelvin for the calculation
    real temp_k = temp + 273.15;
    real Tmin_k = Tmin + 273.15;
    real Tmax_k = Tmax + 273.15;
    real Topt_k = Topt + 273.15;

    // Calculate the numerator of the LRF equation
    real numerator = rmax * (temp_k - Tmax_k) * pow(temp_k - Tmin_k, 2);

    // Calculate the denominator of the LRF equation
    real term1 = (Topt_k - Tmin_k);
    real term2 = term1 * (temp_k - Topt_k);
    real term3 = (Topt_k - Tmax_k) * (Topt_k + Tmin_k - 2.0 * temp_k);
    real denominator = term1 * (term2 - term3);

    // Avoid division by zero and ensure the rate is non-negative
    if (fabs(denominator) < 1e-9) {
      return 0.0;
    } else {
      return fmax(0.0, numerator / denominator);
    }
  }

  /**
   * Non-viability function (quadratic).
   */

  real nonviability(real temp, real a, real Topt, real k) {
    return a * pow(temp - Topt, 2) + k;
  }

  /**
   * Fecundity function.
   */
  real fecundity(real temp, real Topt, real fecmax, real Tb, real K) {
    return pow(20.0, -pow(fabs((temp - Topt) / (Tb / 2.0)), K + 1.0)) * fecmax;
  }
"
stan_funs <- stanvar(scode = stan_funs, block = "functions")
}
