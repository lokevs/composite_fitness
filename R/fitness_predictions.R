#### Preparations ----

# Read thermal performance data
dat <- read.delim("./data/thermal_performance.txt")

# Read oviposition rate/LRS conversion data
ovi_conversion <- read.delim("./data/ovi_rate_conversion.txt")
ovi_conversion$ratio <- ovi_conversion$eggs.1h / ovi_conversion$LRS
ovi_conversion$log.ratio <- log(ovi_conversion$ratio)
ovi_conversion$color <- ifelse(ovi_conversion$origin == "yem", "#ff7f00", 
                               ifelse(ovi_conversion$origin == "bra", "#33a02c", "#1f78b4"))

# Color variable
dat$color <- ifelse(dat$origin == "yem", "#ff7f00", 
                               ifelse(dat$origin == "bra", "#33a02c", "#1f78b4"))

# Read other data
temps <- read.delim("./data/era5_thermal_regimes.txt")

# Packages
library(lubridate)

# Load functions
source("./R/functions.R")

##########

pdf("./output/Fig_S9_LRS_to_ovi_conversion.pdf", width = 2.5, height = 3)
#### Plot the temperature-dependence of LRS to oviposition rate ratio ----

plot(ovi_conversion$eggs.1h / ovi_conversion$LRS ~ ovi_conversion$temperature, 
     pch = 21, 
     bg = ovi_conversion$color,
     cex = 1.5,
     xlim = c(20, 37), ylim = c(0.01, 0.055),
     xlab = "Temperature",
     ylab = "1h oviposition rate / LRS")
curve(fecundity_to_oviposition_rate(fecundity = 1, temperature = x), 20, 37, add = T, lwd = 2)

##########
dev.off()

pdf("./output/Fig_S8_body_size_to_max_fec_conversion.pdf", width = 2.5, height = 3)
#### Plot body size to maximal fecundity  ----

# 29 degrees has highest LRS
aggregate(adult.offspring ~ temperature, data = dat, mean)

# Plot mass to fecundity conversion at 29 degrees
plot(adult.offspring ~ average.mass, data = dat[dat$temperature == 29,], 
     xlab = "Mean adult dry mass (g)", 
     ylab = "Surviving offspring",
     bg = color, pch = 21)
curve(mass_to_fecundity(x, "california"), 0.001, 0.00205, col = unique(dat$color[dat$origin == "ca"]), add = T)
curve(mass_to_fecundity(x, "brazil"), 0.001, 0.00205, col = unique(dat$color[dat$origin == "bra"]), add = T)
curve(mass_to_fecundity(x, "yemen"), 0.001, 0.00205, col = unique(dat$color[dat$origin == "yem"]), add = T)


##########
dev.off()

pdf("./output/Fig_2_fitness_components.pdf", width = 2.5, height = 3)
#### Plot the curves of the separate component functions ----

curve(growth_rate(x, origin = "yemen"), from = 0, to = 50, col = "#ff7f00", lwd = 2, ylab = "Growth rate (g/hour)", xlab = "Temperature", ylim = c(0, 0.0000035))
curve(growth_rate(x, origin = "california"), add = TRUE, col = "#1f78b4", lwd = 2)
curve(growth_rate(x, origin = "brazil"), add = TRUE, col = "#33a02c", lwd = 2)

curve(development_rate(x, origin = "yemen"), from = 0, to = 50, col = "#ff7f00", lwd = 2, ylab = "Development rate (1/hour)", xlab = "Temperature", ylim = c(0, 0.0025))
curve(development_rate(x, origin = "california"), add = TRUE, col = "#1f78b4", lwd = 2)
curve(development_rate(x, origin = "brazil"), add = TRUE, col = "#33a02c", lwd = 2)

curve(fertility_rate(x, origin = "yemen"), from = 0, to = 50, col = "#ff7f00", lwd = 2, ylab = "Fertility rate (P(hour))", xlab = "Temperature", ylim = c(0.95, 1))
curve(fertility_rate(x, origin = "california"), add = TRUE, col = "#1f78b4", lwd = 2)
curve(fertility_rate(x, origin = "brazil"), add = TRUE, col = "#33a02c", lwd = 2)

curve(survival_rate(x, origin = "yemen"), from = 0, to = 50, col = "#ff7f00", lwd = 2, ylab = "Survival rate (P(hour))", xlab = "Temperature", ylim = c(0.7, 1))
curve(survival_rate(x, origin = "california"), add = TRUE, col = "#1f78b4", lwd = 2)
curve(survival_rate(x, origin = "brazil"), add = TRUE, col = "#33a02c", lwd = 2)

##########
dev.off()

pdf("./output/Fig_2_composite_fitness.pdf", width = 3.5, height = 3.5)
#### Plot composite fitness functions over data ----

plot(adult.offspring*dev.rate ~ jitter(temperature),
     data = dat[dat$origin == "yem",],
     bg = color,
     pch = 21,
     lwd = 0.5,
     xlim = c(13, 40), ylim = c(0, 0.22))

points(adult.offspring*dev.rate ~ jitter(temperature),
     data = dat[dat$origin == "bra",],
     bg = color,
     lwd = 0.5,
     pch = 21)
points(adult.offspring*dev.rate ~ jitter(temperature),
     data = dat[dat$origin == "ca",],
     bg = color,
     lwd = 0.5,
     pch = 21)
curve(composite_fitness_function(temperature = x, "yemen", IGR = F), 10, 40, add = T, 
      col = "#ff7f00", lwd = 2)
curve(composite_fitness_function(temperature = x, "brazil", IGR = F), 10, 40, add = T, 
      col = "#33a02c", lwd = 2)
curve(composite_fitness_function(temperature = x, "california", IGR = F), 10, 40, add = T, 
      col = "#1f78b4", lwd = 2)


##########
dev.off()

pdf("./output/Fig_2_temperature_time_series.pdf", width = 3, height = 2.5)
#### Plot temperature time series ----

for(loc in c("california", "brazil", "yemen")){
  
  plot(NA, xlim = c(1, 365*24), ylim = c(-17, 47))
  abline(h = 0, lty = 2)
  if(loc == "california"){ color <- "#1f78b430"}
  else if(loc == "brazil"){ color <- "#33a02c30"}
  else { color <- "#ff7f0030"}
  
  for(site in unique(temps$unique.ID[temps$location == loc])[11:20]){
    
    temp_temps <- subset(temps, unique.ID == site & middle.year == 1)
    lines(y = temp_temps$temperature, x = 1:(365*24), col = color)
  }
  }

##########
dev.off()

pdf("./output/Fig_2_average_temperature_time_series.pdf", width = 3, height = 2.5)
#### Plot average temperature time series ----

plot(NA, xlim = c(1, 365), ylim = c(0, 40))
  abline(h = 0, lty = 2)

for(loc in c("california", "brazil", "yemen")){

    if(loc == "california"){color <- "#1f78b4"; bg_color <- "#1f78b430"}
  else if(loc == "brazil"){color <- "#33a02c"; bg_color <- "#33a02c30"}
  else {color <- "#ff7f00"; bg_color <- "#ff7f0030"}

  temp_temps <- subset(temps, location == loc & middle.year == 1)
  daily_means <- aggregate(temperature ~ format(as_datetime(date.time), "%m-%d"), temp_temps, "mean")
  daily_variances <- aggregate(temperature ~ format(as_datetime(date.time), "%m-%d") + unique.ID, temp_temps, "sd")
  colnames(daily_variances) <- c("dates", "unique.ID", "temperature")
  daily_mean_variances <- aggregate(temperature ~ dates, daily_variances, "mean")
  
  polygon(x = c(1:365, 365:1),
    y = c(daily_means$temperature+daily_mean_variances$temperature, rev(daily_means$temperature-(daily_mean_variances$temperature))),
    col = bg_color,
    border = F)
  
  lines(y = daily_means$temperature, x = 1:365, col = color, lwd = 1.5)

}

##########
dev.off()

#### Predict fitness in variable settings ----

# Define the warming scenarios (mean temperature increase)
warming_scenarios <- seq(from = 0, to = 5.5, by = 0.1)

# Get the unique geographic locations from the temperature data
unique_locations <- unique(temps$unique.ID)

# Initialize an empty list to store the results from each iteration.
# Using a list is more memory-efficient than repeatedly binding rows to a data frame.
all_results_list <- list()

# Run simulations
for (location_id in unique_locations) {
  
  # Subset the temperature data for the current location
  location_data <- subset(temps, unique.ID == location_id)
  
  # Ensure there is at least a three years of hourly data (365 * 24 = 8760 hours)
  if (nrow(location_data) < 365 * 24 * 3 - 24) {
    warning(paste("Skipping location", location_id, "- insufficient data (< 3 years)."))
    next # Skip to the next location
  }
  
  # Extract base temperature data and coordinates
  base_temperatures <- location_data$temperature
  lat <- location_data$latitude[1]
  lon <- location_data$longitude[1]
  location_name <- location_data$location[1] 
  

    # Loop through each warming scenario
    for (warming_level in warming_scenarios) {
      
      # Apply the warming shift to the base temperatures for the current scenario
      warmed_temperatures <- base_temperatures + warming_level
      
      # Run simulation and fitness calculation
      summary_output <- calculate_fitness(
        temperatures = warmed_temperatures,
        origin = location_name, 
        plot = FALSE
      )
      
      # Add a row for the "separate" (full simulation) fitness result
      all_results_list[[length(all_results_list) + 1]] <- data.frame(
        location = location_name,
        origin = location_name,
        latitude = lat,
        longitude = lon,
        unique.ID = location_id,
        mean.fitness = summary_output$lambda.separate.ari ,
        geomean.fitness = summary_output$lambda.separate.geo ,
        function.type = "separate",
        warming = warming_level
      )
      
      # Add a row for the "composite" (direct function) fitness result
      all_results_list[[length(all_results_list) + 1]] <- data.frame(
        location = location_name,
        origin = location_name,
        latitude = lat,
        longitude = lon,
        unique.ID = location_id,
        mean.fitness = summary_output$lambda.composite.ari,
        geomean.fitness = summary_output$lambda.composite.geo,
        function.type = "composite",
        warming = warming_level
      )
      
    } 
  
  # Provide progress update in the console
  cat(paste("  - Completed processing for location:", location_id, "\n"))
  
} 


# Combine the list of results into a data frame
all_fitness_simus <- do.call(rbind, all_results_list)

write.table(all_fitness_simus, "./data/simulated_fitness.txt", sep = "\t", row.names = F, quote = F)

##########

#### Run simulation under controlled yearly/daily fluctuation amplitudes ####
generate_temp_series <- function(mean_temp, yearly_amp, daily_amp) {
  
  # Define vectors
  total_hours <- 3 * 365 * 24
  t_hours <- 0:(total_hours - 1) 
  
  # Yearly component
  period_yearly_hours <- 365 * 24
  yearly_component <- yearly_amp * (-cos(2 * pi * t_hours / period_yearly_hours))
  
  # Daily component
  period_daily_hours <- 24
  daily_component <- daily_amp * (-cos(2 * pi * t_hours / period_daily_hours))
  
  # Combine
  # The final temperature is the sum of the mean and the two components
  temperatures <- mean_temp + yearly_component + daily_component
  
  return(temperatures)
}


# Set the range of parameters 
mean_temps <- seq(from = 10, to = 30, by = 5)  
yearly_amps <- seq(from = 0, to = 15, by = 0.5) 
daily_amps <- seq(from = 0, to = 15, by = 0.5)    
origins_to_run <- c("california", "brazil", "yemen")

# For counter
total_iterations <- length(mean_temps) * length(yearly_amps) * length(daily_amps) * length(origins_to_run)

# Initialize list
all_results_list <- list()
iter <- 1

# Loop
for (mean_val in mean_temps) {
  for (yearly_val in yearly_amps) {
    for (daily_val in daily_amps) {
      for (origin in origins_to_run) {
        
        # 1. Generate Temperatures
        base_temperatures <- generate_temp_series(mean_val, yearly_val, daily_val)
        
        # Main function
        summary_output <- calculate_fitness(
          temperatures = base_temperatures,
          origin = origin, 
          plot = FALSE
        )
        
        # 2. Run Models Directly
        sep_res <- separate_fitness_function(base_temperatures, origin)
        
        # Slice Year 2 for composite
        year_hours <- 365 * 24
        temps_year2 <- if(length(base_temperatures) >= (2 * year_hours)) {
          base_temperatures[(year_hours + 1):(2 * year_hours)] 
        } else { 
          base_temperatures 
        }
        comp_res <- composite_fitness_function(temps_year2, origin, return_components = TRUE)
        
        # 3. Calculate Totals (Reconstruct from Daily ^ DevTimeDays)
        # S_total = (S_daily)^(Dev_Days)
        S_tot_sep_vec <- sep_res$survival.daily ^ sep_res$dev.time
        S_tot_comp_vec <- comp_res$survival.daily ^ comp_res$dev.time
        
        # F_total = (F_daily)^(Dev_Days)
        F_tot_sep_vec <- sep_res$fertility.daily ^ sep_res$dev.time
        F_tot_comp_vec <- comp_res$fertility.daily ^ comp_res$dev.time
        
        # 4. Save Results
        all_results_list[[length(all_results_list) + 1]] <- data.frame(
          origin = origin,
          mean_temp = mean_val,
          yearly_amplitude = yearly_val,
          daily_amplitude = daily_val,
          
          # Traits: lambda
          mean_fitness_separate_ari = summary_output$lambda.separate.ari,
          mean_fitness_composite_ari = summary_output$lambda.composite.ari,
          mean_fitness_separate_geo = summary_output$lambda.separate.geo,
          mean_fitness_composite_geo = summary_output$lambda.composite.geo,
          
          # Traits: R0
          R0_sep = mean(sep_res$R0, na.rm = TRUE),
          R0_comp = mean(comp_res$R0, na.rm = TRUE),
          
          # Traits: Components for Gen Time
          dev_time_sep = mean(sep_res$dev.time, na.rm = TRUE),
          dev_time_comp = mean(comp_res$dev.time, na.rm = TRUE),
          ovi_time_sep = mean(sep_res$ovi.time.days, na.rm = TRUE),
          ovi_time_comp = mean(comp_res$ovi.time.days, na.rm = TRUE),
          
          # Traits: Mass
          mass_sep = mean(sep_res$mass, na.rm = TRUE),
          mass_comp = mean(comp_res$mass, na.rm = TRUE),
          
          # Traits: TOTAL Survival
          surv_tot_sep = mean(S_tot_sep_vec, na.rm = TRUE),
          surv_tot_comp = mean(S_tot_comp_vec, na.rm = TRUE),
          
          # Traits: TOTAL Fertility
          fert_tot_sep = mean(F_tot_sep_vec, na.rm = TRUE),
          fert_tot_comp = mean(F_tot_comp_vec, na.rm = TRUE)
        )
        
        if(iter %% 100 == 0) cat(paste("- ", iter, " / ", total_iterations, "\n"))
        iter = iter + 1
      }
    } 
  } 
}

# Combine 
all_fitness_simus <- do.call(rbind, all_results_list)

# Sum 
all_fitness_simus$fitness.diff.ari <- ifelse(all_fitness_simus$mean_fitness_composite_ari + all_fitness_simus$mean_fitness_separate_ari == 0, 0, 
                                             (all_fitness_simus$mean_fitness_composite_ari - all_fitness_simus$mean_fitness_separate_ari) / (all_fitness_simus$mean_fitness_composite_ari + all_fitness_simus$mean_fitness_separate_ari))


# Product (scaled relative difference between -1 and 1, tanh method to avoid overly large numbers)
log_diff <- log(all_fitness_simus$mean_fitness_composite_geo) - log(all_fitness_simus$mean_fitness_separate_geo)
val_input <- 0.5 * 365 * 24 * log_diff
all_fitness_simus$fitness.diff.geo <- tanh(val_input)
zero_mask <- all_fitness_simus$mean_fitness_composite_geo == 0 & all_fitness_simus$mean_fitness_separate_geo == 0
all_fitness_simus$fitness.diff.geo[zero_mask] <- 0

write.table(all_fitness_simus, "./data/simulated_fitness_controlled_settings.txt", sep = "\t", row.names = F, quote = F)

##########

pdf("./output/Fig_S4_fitness_diff_heatmaps.pdf", width = 2.5, height = 3)
#### Plot heatmaps of the controlled setting fitness diffs rate summation ----

simulated_fitness_controlled <- read.delim("./data/simulated_fitness_controlled_settings.txt")

for(mean_temps in unique(simulated_fitness_controlled$mean_temp)){
  
  # Create color variable
  data_subset<- aggregate(fitness.diff.ari ~ yearly_amplitude + daily_amplitude + mean_temp, simulated_fitness_controlled, "mean")
  data_subset$color.map <- round((data_subset[,4]* 50 + 51))
  
  # Subset for temp 
  data_subset <- data_subset[data_subset$mean_temp == mean_temps, ]
  
  # Create matrix with colors
  heatmap_matrix <- round(tapply(data_subset$color.map, 
                                 list(data_subset$daily_amplitude, data_subset$yearly_amplitude), 
                                 FUN = mean))
  heatmap_matrix[] <- colorRampPalette(c("darkcyan", "white", "red"))(101)[heatmap_matrix]
  
  # Set up empty plot
  plot(NULL, 
       xlim = c(1, ncol(heatmap_matrix)), 
       ylim = c(1, nrow(heatmap_matrix)),
       xlab = "Yearly amplitude (°C)", 
       ylab = "Daily amplitude (°C)",
       xaxt = "n",
       yaxt = "n",
       main = paste0("Scaled fit. diff, ", mean_temps))
  
  # Add axes 
  axis(1, at = 1:ncol(heatmap_matrix), labels = colnames(heatmap_matrix))
  axis(2, at = 1:nrow(heatmap_matrix), labels = rownames(heatmap_matrix), las = 1)
  
  # Plot the matrix of colors
  rasterImage(heatmap_matrix[nrow(heatmap_matrix):1, ], 
              xleft = 0.5, 
              ytop = nrow(heatmap_matrix) + 0.5, 
              xright = ncol(heatmap_matrix) + 0.5, 
              ybottom = 0.5,
              interpolate = F)
}

##########

#### Plot heatmaps of the controlled setting fitness diffs geometric mean ----

simulated_fitness_controlled <- read.delim("./data/simulated_fitness_controlled_settings.txt")

for(mean_temps in unique(simulated_fitness_controlled$mean_temp)){
  
  # Create color variable
  data_subset<- aggregate(fitness.diff.geo ~ yearly_amplitude + daily_amplitude + mean_temp, simulated_fitness_controlled, "mean")
  data_subset$color.map <- round((data_subset[,4]* 50 + 51))
  
  # Subset for temp 
  data_subset <- data_subset[data_subset$mean_temp == mean_temps, ]
  
  # Create matrix with colors
  heatmap_matrix <- round(tapply(data_subset$color.map, 
                                 list(data_subset$daily_amplitude, data_subset$yearly_amplitude), 
                                 FUN = mean))
  heatmap_matrix[] <- colorRampPalette(c("darkcyan", "white", "red"))(101)[heatmap_matrix]
  
  # Set up empty plot
  plot(NULL, 
       xlim = c(1, ncol(heatmap_matrix)), 
       ylim = c(1, nrow(heatmap_matrix)),
       xlab = "Yearly amplitude (°C)", 
       ylab = "Daily amplitude (°C)",
       xaxt = "n",
       yaxt = "n",
       main = paste0("Scaled fit. diff, ", mean_temps))
  
  # Add axes 
  axis(1, at = 1:ncol(heatmap_matrix), labels = colnames(heatmap_matrix))
  axis(2, at = 1:nrow(heatmap_matrix), labels = rownames(heatmap_matrix), las = 1)
  
  # Plot the matrix of colors
  rasterImage(heatmap_matrix[nrow(heatmap_matrix):1, ], 
              xleft = 0.5, 
              ytop = nrow(heatmap_matrix) + 0.5, 
              xright = ncol(heatmap_matrix) + 0.5, 
              ybottom = 0.5,
              interpolate = F)
}

##########
dev.off()

pdf("./output/Fig_5_trait_heatmaps.pdf", width = 2.5, height = 3)
#### Plot Trait Differences ----

# 1. Read data
sim_data <- read.delim("./data/simulated_fitness_controlled_settings.txt")

# 2. Subset temperatures
target_temps <- c(10, 15, 20, 25, 30)
sim_data <- sim_data[sim_data$mean_temp %in% target_temps, ]

# 3. Calculate Arithmetic Differences (Red/Positive = "Better" for Composite)

# A. Gen Time: Shorter is better.
# Positive difference means Sep > Comp (Comp is shorter/better)
gen_sep <- sim_data$dev_time_sep + (sim_data$ovi_time_sep / 2)
gen_comp <- sim_data$dev_time_comp + (sim_data$ovi_time_comp / 2)
sim_data$diff_gen <- (gen_sep - gen_comp)

# B. Max Fecundity: Higher is better. (Replacing Mass)
# We must calculate Max Fecundity from Mass + Origin manually here to vectorize it
# Intercepts from functions.R: CA=6.9582461, BRA=7.1822847, YEM=7.2572643
intercepts <- c("california" = 6.9582461, "brazil" = 7.1822847, "yemen" = 7.2572643)
int_vec <- intercepts[sim_data$origin]

# Formula: exp(intercept + 0.4440254 * log(mass))
fec_max_sep <- exp(int_vec + 0.4440254 * log(sim_data$mass_sep))
fec_max_comp <- exp(int_vec + 0.4440254 * log(sim_data$mass_comp))

# Difference: Comp - Sep (Positive means Comp is higher/better)
sim_data$diff_fec_max <- (fec_max_comp - fec_max_sep)

# C. Survival (TOTAL): Higher is better.
sim_data$diff_surv <- (sim_data$surv_tot_comp - sim_data$surv_tot_sep)

# D. Fertility (TOTAL): Higher is better.
sim_data$diff_fert <- (sim_data$fert_tot_comp - sim_data$fert_tot_sep)

# Handle NaN cases
cols_to_fix <- c("diff_gen", "diff_fec_max", "diff_surv", "diff_fert")
for(c in cols_to_fix) {
  if(c %in% colnames(sim_data)) {
    sim_data[[c]][is.nan(sim_data[[c]])] <- 0
  }
}

# Define traits
traits_to_plot <- list(
  list(col = "diff_gen", title = "Gen. Time (Days Diff)"),
  list(col = "diff_fec_max", title = "Max Fecundity (Eggs Diff)"),
  list(col = "diff_surv", title = "Total Survival (Prob Diff)"),
  list(col = "diff_fert", title = "Total Fertility (Prob Diff)")
)

# Define divergent palette
n_colors <- 101
col_palette <- colorRampPalette(c("darkcyan", "white", "red"))(n_colors)

# Special colors for Infinity
col_pos_inf <- "#8B0000" # Dark Red
col_neg_inf <- "#008B8B" # Dark Cyan

for (trait in traits_to_plot) {
  
  # --- DYNAMIC SCALE (Per Row) ---
  global_max_abs <- 0
  for(t in target_temps){
    sub_dat <- sim_data[sim_data$mean_temp == t, ]
    vals <- sub_dat[[trait$col]]
    finite_vals <- vals[is.finite(vals)]
    if(length(finite_vals) > 0){
      # Aggregation to match pixel resolution
      agg <- aggregate(finite_vals ~ sub_dat$yearly_amplitude[is.finite(vals)] + sub_dat$daily_amplitude[is.finite(vals)], FUN=mean)
      global_max_abs <- max(global_max_abs, max(abs(agg[,3]), na.rm = TRUE))
    }
  }
  
  if(global_max_abs == 0) global_max_abs <- 1
  z_lims <- c(-global_max_abs, global_max_abs)
  # -------------------------------
  
  for (temp in target_temps) {
    
    # Subset
    data_subset <- sim_data[sim_data$mean_temp == temp, ]
    
    # Aggregate
    agg_data <- aggregate(data_subset[[trait$col]] ~ yearly_amplitude + daily_amplitude, data_subset, function(x) {
      if(any(is.infinite(x))) {
        return(x[is.infinite(x)][1]) 
      } else {
        return(mean(x, na.rm=TRUE))
      }
    })
    colnames(agg_data)[3] <- "value"
    
    # Matrix
    heatmap_matrix <- tapply(agg_data$value, 
                             list(agg_data$daily_amplitude, agg_data$yearly_amplitude), 
                             mean)
    
    # Normalize
    normalized_vals <- (heatmap_matrix - z_lims[1]) / (z_lims[2] - z_lims[1])
    
    # Clamp
    normalized_vals[is.finite(normalized_vals) & normalized_vals < 0] <- 0
    normalized_vals[is.finite(normalized_vals) & normalized_vals > 1] <- 1
    normalized_vals[is.nan(normalized_vals)] <- 0.5 
    
    # Colors
    color_idx <- round(normalized_vals * (n_colors - 1)) + 1
    base_colors <- col_palette[color_idx]
    
    final_colors <- base_colors
    final_colors[which(heatmap_matrix == Inf)] <- col_pos_inf
    final_colors[which(heatmap_matrix == -Inf)] <- col_neg_inf
    
    # Reshape
    col_mat <- matrix(final_colors, nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix))
    col_mat_flipped <- col_mat[nrow(col_mat):1, ]
    
    # Plot
    plot(NULL, 
         xlim = c(1, ncol(heatmap_matrix)), 
         ylim = c(1, nrow(heatmap_matrix)),
         xaxt = "n", yaxt = "n", bty = "n",
         main = "", xlab = "", ylab = "")
    
    if (trait$col == traits_to_plot[[1]]$col) {
      mtext(paste0(temp, "°C"), side = 3, line = 0.5, cex = 0.8)
    }
    
    rasterImage(col_mat_flipped, 
                xleft = 0.5, ytop = nrow(heatmap_matrix) + 0.5, 
                xright = ncol(heatmap_matrix) + 0.5, ybottom = 0.5,
                interpolate = FALSE)
    
    if (temp == target_temps[1]) {
      axis(2, at = 1:nrow(heatmap_matrix), labels = rownames(heatmap_matrix), las = 1, cex.axis = 0.7)
      mtext(trait$title, side = 2, line = 2.5, cex = 0.6, font = 2)
    }
    
    if (trait$col == traits_to_plot[[length(traits_to_plot)]]$col) {
      axis(1, at = 1:ncol(heatmap_matrix), labels = colnames(heatmap_matrix), cex.axis = 0.7)
    }
    
    box()
  }
}

mtext("Yearly Amplitude", side = 1, outer = TRUE, line = 2)
mtext("Daily Amplitude", side = 2, outer = TRUE, line = 2)

##########

#### Plot Custom Gradient Bars for Each Trait ----

sim_data <- read.delim("./data/simulated_fitness_controlled_settings.txt")


gen_sep <- sim_data$dev_time_sep + (sim_data$ovi_time_sep / 2)
gen_comp <- sim_data$dev_time_comp + (sim_data$ovi_time_comp / 2)
sim_data$diff_gen <- (gen_sep - gen_comp)

intercepts <- c("california" = 6.9582461, "brazil" = 7.1822847, "yemen" = 7.2572643)
int_vec <- intercepts[sim_data$origin]
fec_max_sep <- exp(int_vec + 0.4440254 * log(sim_data$mass_sep))
fec_max_comp <- exp(int_vec + 0.4440254 * log(sim_data$mass_comp))
sim_data$diff_fec_max <- (fec_max_comp - fec_max_sep)

sim_data$diff_surv <- (sim_data$surv_tot_comp - sim_data$surv_tot_sep)
sim_data$diff_fert <- (sim_data$fert_tot_comp - sim_data$fert_tot_sep)

traits_to_plot <- list(
  list(col = "diff_gen", title = "Gen. Time (Days Diff)"),
  list(col = "diff_fec_max", title = "Max Fecundity (Eggs Diff)"),
  list(col = "diff_surv", title = "Total Survival (Prob Diff)"),
  list(col = "diff_fert", title = "Total Fertility (Prob Diff)")
)

legend_colors <- colorRampPalette(c("darkcyan", "white", "red"))(101)
color_matrix <- as.raster(matrix(legend_colors, nrow = 1))

for (trait in traits_to_plot) {
  
  actual_min <- Inf
  actual_max <- -Inf
  global_max_abs <- 0
  
  for(t in target_temps){
    sub_dat <- sim_data[sim_data$mean_temp == t, ]
    vals <- sub_dat[[trait$col]]
    finite_vals <- vals[is.finite(vals)]
    if(length(finite_vals) > 0){
      agg <- aggregate(finite_vals ~ sub_dat$yearly_amplitude[is.finite(vals)] + sub_dat$daily_amplitude[is.finite(vals)], FUN=mean)
      actual_min <- min(actual_min, min(agg[,3], na.rm = TRUE))
      actual_max <- max(actual_max, max(agg[,3], na.rm = TRUE))
      global_max_abs <- max(global_max_abs, max(abs(agg[,3]), na.rm = TRUE))
    }
  }
  if(global_max_abs == 0) global_max_abs <- 1
  if(actual_min == Inf) actual_min <- -1
  if(actual_max == -Inf) actual_max <- 1
  
  # Empty plot
  plot(NA, xlim = c(actual_min, actual_max), ylim = c(0, 2), 
       xaxt = 'n', yaxt = 'n', bty = 'n', xlab = "", ylab = "")
  
  # Plot gradient raster
  rasterImage(color_matrix, 
              xleft = -global_max_abs, 
              xright = global_max_abs, 
              ybottom = 0.5, 
              ytop = 1)
  
  tick_vals <- c(actual_min, actual_max)
  if(actual_min < 0 && actual_max > 0) tick_vals <- c(actual_min, 0, actual_max)
  
  # Tick labels
  text(x = tick_vals, 
       y = 0.4,       
       pos = 1,        
       labels = round(tick_vals, 2),
       cex = 1.2)      
  
  # Title
  text(x = mean(c(actual_min, actual_max)), y = 1.2, pos = 3, 
       labels = trait$title,
       cex = 1.2,
       font = 2) 
}

##########
dev.off()
