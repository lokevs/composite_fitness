#### Preparations ----

# Read data
dat <- read.delim("./data/simulated_fitness.txt")
temps <- read.delim("./data/era5_thermal_regimes.txt")
source("./R/functions.R")

# Make it wide
dat_wide <- reshape(dat,
                    idvar = c("location", "origin", "latitude", "longitude", "unique.ID", "warming"),
                    timevar = "function.type",
                    v.names = c("mean.fitness", "geomean.fitness"),
                    direction = "wide")

# Colors and the effects of warming
dat_wide$color <- ifelse(dat_wide$location == "california", "#1f78b4", ifelse(dat_wide$location == "yemen", "#ff7f00", "#33a02c"))

# Calculate fitness differences

# Sum (scaled relative difference between -1 and 1)
dat_wide$fitness.diff.ari <- ifelse(dat_wide$mean.fitness.composite + dat_wide$mean.fitness.separate == 0, 0, 
                                    (dat_wide$mean.fitness.composite - dat_wide$mean.fitness.separate) / (dat_wide$mean.fitness.composite + dat_wide$mean.fitness.separate))

# Product (scaled relative difference between -1 and 1, tanh method to avoid overly large numbers)
log_diff <- log(dat_wide$geomean.fitness.composite) - log(dat_wide$geomean.fitness.separate)
val_input <- 0.5 * 365 * 24 * log_diff
dat_wide$fitness.diff.geo <- tanh(val_input)
zero_mask <- dat_wide$geomean.fitness.composite == 0 & dat_wide$geomean.fitness.separate == 0
dat_wide$fitness.diff.geo[zero_mask] <- 0


# Calculate some metrics
# Sort to ensure time-series operations work correctly (Critical for diff/cummax)
dat_wide <- dat_wide[order(dat_wide$unique.ID, dat_wide$warming), ]

# Helper functions for ave()
# Negative warming flag
# Returns 1 if fitness is declining OR (fitness is 0 but was previously > 0)
calc_neg_warming <- function(x) {
  change <- c(0, diff(x))
  is_declining <- change < 0
  is_crashed   <- (x == 0) & (cummax(x) > 0)
  as.integer(is_declining | is_crashed)
}

# Apply metrics for 'separate (ari)
dat_wide$negative.warming.separate.ari     <- ave(dat_wide$mean.fitness.separate, dat_wide$unique.ID, FUN = calc_neg_warming)

# Apply metrics for composite (ari)
dat_wide$negative.warming.composite.ari    <- ave(dat_wide$mean.fitness.composite, dat_wide$unique.ID, FUN = calc_neg_warming)

# Apply metrics for separate (geo)
dat_wide$negative.warming.separate.geo     <- ave(dat_wide$geomean.fitness.separate, dat_wide$unique.ID, FUN = calc_neg_warming)

# Apply metrics for composite (geo)
dat_wide$negative.warming.composite.geo    <- ave(dat_wide$geomean.fitness.composite, dat_wide$unique.ID, FUN = calc_neg_warming)

# Habitable flags 
dat_wide$habitable.separate.ari <- as.integer(dat_wide$mean.fitness.separate >= 1)
dat_wide$habitable.composite.ari <- as.integer(dat_wide$mean.fitness.composite >= 1)
dat_wide$habitable.separate.geo <- as.integer(dat_wide$geomean.fitness.separate >= 1)
dat_wide$habitable.composite.geo <- as.integer(dat_wide$geomean.fitness.composite >= 1)

# Extinction flags 
dat_wide$extinction.separate.ari <- ifelse(dat_wide$mean.fitness.separate == 0, 1, 0)
dat_wide$extinction.composite.ari <- ifelse(dat_wide$mean.fitness.composite == 0, 1, 0)
dat_wide$extinction.separate.geo <- ifelse(dat_wide$geomean.fitness.separate == 0, 1, 0)
dat_wide$extinction.composite.geo <- ifelse(dat_wide$geomean.fitness.composite == 0, 1, 0)

# Libraries
library(maps)

##########

pdf("./output/Fig_S6_site_specific_mean_fitness.pdf", width = 2.5, height = 3)
#### Plot of mean fitness as a function of warming ----

# Helper function for linear interpolation of x where y crosses 1
get_crossing_x <- function(x1, y1, x2, y2, target_y = 1) {
  if (abs(y2 - y1) < 1e-9) return(mean(c(x1, x2)))
  m <- (y2 - y1) / (x2 - x1)
  c_int <- y1 - m * x1
  return((target_y - c_int) / m)
}


# Define plotting parameters to loop through to avoid code repetition
plot_configs <- list(
  list(data_col = "mean.fitness.separate", hab_col = "habitable.separate.ari", title_suffix = "separate functions, arithmetic mean", ylab = "Mean fitness"),
  list(data_col = "mean.fitness.composite", hab_col = "habitable.composite.ari", title_suffix = "composite functions, arithmetic mean", ylab = "Lambda daily"),
  list(data_col = "geomean.fitness.separate", hab_col = "habitable.separate.geo", title_suffix = "separate functions, geometric mean", ylab = "Mean fitness"),
  list(data_col = "geomean.fitness.composite", hab_col = "habitable.composite.geo", title_suffix = "composite functions, geometric mean", ylab = "Lambda daily")
)

for (cfg in plot_configs) {
  for (loc in c("california", "brazil", "yemen")) {
    
    plot(NA, ylim = c(0, 2.05), xlim = c(0, 5), 
         xlab = "Warming (°C)", ylab = cfg$ylab, 
         main = paste0(loc, ", ", cfg$title_suffix), 
         yaxt = "n")
    
    # Custom Y-axis
    axis(2, c(0, 1, 1.007^(4*24)), c(0, 1, 1.007))
    abline(h = 1, lty = 2, col = "gray60")
    
    # Filter data for location
    loc_data <- dat_wide[dat_wide$location == loc & dat_wide$warming <= 5.1, ]
    
    for (i in unique(loc_data$unique.ID)) {
      
      sub_dat <- loc_data[loc_data$unique.ID == i, ]
      sub_dat <- sub_dat[order(sub_dat$warming), ]
      
      # Extract vectors for easier handling
      x <- sub_dat$warming
      # Raise to power as per your original logic
      y <- sub_dat[[cfg$data_col]]^(4*24) 
      habitable <- sub_dat[[cfg$hab_col]]
      sub_col <- unique(sub_dat$color)
      
      # --- 1. Draw Lines (Segment by Segment) ---
      if (length(x) > 1) {
        for (j in 1:(length(x) - 1)) {
          
          # Check if declining
          is_declining <- y[j+1] < y[j]
          
          # Determine Style based ONLY on decline status
          if (is_declining) {
            # Declining: Full color, thick line
            use_col <- sub_col
            use_lwd <- 1 
          } else {
            # Increasing or flat: Faint color, thin line
            use_col <- paste0(sub_col, "35") 
            use_lwd <- 0.5
          }
          
          lines(x[j:(j+1)], y[j:(j+1)], col = use_col, lwd = use_lwd)
        }
      }
      
      # --- 2. Draw Transition Points (Interpolated) ---
      # Transitions happen when habitable status changes
      trans <- diff(habitable)
      recov_indices <- which(trans == 1) # 0 -> 1
      loss_indices <- which(trans == -1) # 1 -> 0
      
      # Plot Recovery Points (0 -> 1)
      if (length(recov_indices) > 0) {
        for (idx in recov_indices) {
          # Interpolate x where y crosses 1 between idx and idx+1
          x_cross <- get_crossing_x(x[idx], y[idx], x[idx+1], y[idx+1], target_y = 1)
          # Ensure x_cross is within bounds (handles potential numerical issues)
          x_cross <- max(x[idx], min(x_cross, x[idx+1]))
          points(x_cross, 1, pch = 21, col = sub_col, bg = "white", cex = 1.2)
        }
      }
      
      # Plot Loss Points (1 -> 0)
      if (length(loss_indices) > 0) {
        for (idx in loss_indices) {
          x_cross <- get_crossing_x(x[idx], y[idx], x[idx+1], y[idx+1], target_y = 1)
          x_cross <- max(x[idx], min(x_cross, x[idx+1]))
          points(x_cross, 1, pch = 4, col = sub_col, lwd = 1.3, cex = 1.2)
        }
      }
    }
  }
}


# Plot number of sites in decline
aggregated_separate <- aggregate(negative.warming.separate.ari ~ origin + color + warming, data = dat_wide, "mean")
aggregated_composite <- aggregate(negative.warming.composite.ari ~ origin + color + warming, data = dat_wide, "mean")

for(origin in unique(aggregated_separate$origin)){
  
  data_separate <- aggregated_separate[aggregated_separate$warming %in% c(0, 2.5, 5) & aggregated_separate$origin == origin,]
  data_composite <- aggregated_composite[aggregated_composite$warming %in% c(0, 2.5, 5) & aggregated_composite$origin == origin,]
  
  plot_matrix <- rbind(data_separate$negative.warming.separate.ari, 
                       data_composite$negative.warming.composite.ari)
  print(plot_matrix)
  print(plot_matrix[2,] / plot_matrix[1,]-1)
  
  # Create the bar chart
  barplot(plot_matrix,
          beside = TRUE,
          main = origin,
          xlab = "Warming (°C)",
          ylab = "Proportion of sites where additional warming is detrimental",
          ylim = c(0, 1),
          col = unique(data_separate$color),
          density = c(NA, 17),
          bty = "o",
          legend.text = c("Separate", "Composite"))
}  


aggregated_separate <- aggregate(negative.warming.separate.geo ~ origin + color + warming, data = dat_wide, "mean")
aggregated_composite <- aggregate(negative.warming.composite.geo ~ origin + color + warming, data = dat_wide, "mean")

for(origin in unique(aggregated_separate$origin)){
  
  data_separate <- aggregated_separate[aggregated_separate$warming %in% c(0, 2.5, 5) & aggregated_separate$origin == origin,]
  data_composite <- aggregated_composite[aggregated_composite$warming %in% c(0, 2.5, 5) & aggregated_composite$origin == origin,]
  
  plot_matrix <- rbind(data_separate$negative.warming.separate.geo, 
                       data_composite$negative.warming.composite.geo)
  print(plot_matrix)
  print(plot_matrix[2,] / plot_matrix[1,]-1)
  
  # Create the bar chart
  barplot(plot_matrix,
          beside = TRUE,
          main = origin,
          xlab = "Warming (°C)",
          ylab = "Proportion of sites where additional warming is detrimental",
          ylim = c(0, 1),
          col = unique(data_separate$color),
          density = c(NA, 17),
          bty = "o",
          legend.text = c("Separate", "Composite"))
}  


# Plot number of sites going extinct (lambda = 0)
aggregated_separate <- aggregate(extinction.separate.ari ~ origin + color + warming, data = dat_wide, "mean")
aggregated_composite <- aggregate(extinction.composite.ari ~ origin + color + warming, data = dat_wide, "mean")

for(origin in unique(aggregated_separate$origin)){
  
  data_separate <- aggregated_separate[aggregated_separate$warming %in% c(0, 2.5, 5) & aggregated_separate$origin == origin,]
  data_composite <- aggregated_composite[aggregated_composite$warming %in% c(0, 2.5, 5) & aggregated_composite$origin == origin,]
  
  plot_matrix <- rbind(data_separate$extinction.separate.ari, 
                       data_composite$extinction.composite.ari)
  print(plot_matrix)
  print(plot_matrix[2,] / plot_matrix[1,]-1)
  
  # Create the bar chart
  barplot(plot_matrix,
          beside = TRUE,
          main = origin,
          xlab = "Warming (°C)",
          ylab = "Proportion of sites where additional warming is detrimental",
          ylim = c(0, 1),
          col = unique(data_separate$color),
          density = c(NA, 17),
          bty = "o",
          legend.text = c("Separate", "Composite"))
}  


aggregated_separate <- aggregate(extinction.separate.geo ~ origin + color + warming, data = dat_wide, "mean")
aggregated_composite <- aggregate(extinction.composite.geo ~ origin + color + warming, data = dat_wide, "mean")

for(origin in unique(aggregated_separate$origin)){
  
  data_separate <- aggregated_separate[aggregated_separate$warming %in% c(0, 2.5, 5) & aggregated_separate$origin == origin,]
  data_composite <- aggregated_composite[aggregated_composite$warming %in% c(0, 2.5, 5) & aggregated_composite$origin == origin,]
  
  plot_matrix <- rbind(data_separate$extinction.separate.geo, 
                       data_composite$extinction.composite.geo)
  print(plot_matrix)
  print(plot_matrix[2,] / plot_matrix[1,]-1)
  
  # Create the bar chart
  barplot(plot_matrix,
          beside = TRUE,
          main = origin,
          xlab = "Warming (°C)",
          ylab = "Proportion of sites where additional warming is detrimental",
          ylim = c(0, 1),
          col = unique(data_separate$color),
          density = c(NA, 17),
          bty = "o",
          legend.text = c("Separate", "Composite"))
}  

##########
dev.off()

pdf("./output/Fig_3_site_specific_relative_fitness_diff.pdf", width = 2.5, height = 3)
#### Plot of relative fitness differences ----

# Plot arithmetic mean fitness 
for(loc in c("california", "brazil", "yemen")){
  plot(NA, ylim = c(-1, 1), xlim = c(0, 5), xlab = "Warming (°C)", ylab = "Scaled difference in geometric fitness", main = paste0(loc))
  abline(h = 0, lty = 2)
  for(i in unique(dat_wide$unique.ID[dat_wide$location == loc])){
    lines(fitness.diff.ari + runif(1, -0.05, 0.05) ~ warming, 
          data = dat_wide[dat_wide$unique.ID == i & dat_wide$warming <= 5,], 
          col = paste0(color, "35"), lwd = 0.5)
  }
  if(loc == "california") col2 <- "#16547e"
  if(loc == "brazil") col2 <- "#24701f"
  if(loc == "yemen") col2 <- "#b35900"
  
  lines(fitness.diff.ari ~ warming, 
        data = aggregate(fitness.diff.ari ~ warming, dat_wide[dat_wide$location == loc & dat_wide$warming <= 5,], mean), 
        col = col2, lwd = 1.5)
  }

# Plot geometric mean fitness
for(loc in c("california", "brazil", "yemen")){
  plot(NA, ylim = c(-1, 1), xlim = c(0, 5), xlab = "Warming (°C)", ylab = "Scaled difference in geometric fitness", main = paste0(loc))
  abline(h = 0, lty = 2)
  for(i in unique(dat_wide$unique.ID[dat_wide$location == loc])){
    lines(fitness.diff.geo + runif(1, -0.05, 0.05) ~ warming, 
          data = dat_wide[dat_wide$unique.ID == i & dat_wide$warming <= 5,], 
          col = paste0(color, "35"), lwd = 0.5)
  }
  if(loc == "california") col2 <- "#16547e"
  if(loc == "brazil") col2 <- "#24701f"
  if(loc == "yemen") col2 <- "#b35900"

  lines(fitness.diff.geo ~ warming, 
        data = aggregate(fitness.diff.geo ~ warming, dat_wide[dat_wide$location == loc & dat_wide$warming <= 5,], mean), 
        col = col2, lwd = 1.5)
  
}

##########
dev.off()

pdf("./output/Fig_4_example_sites_hourly.pdf", width = 3, height = 2.5)
#### Plot of six example sites with different and their hourly fitness ----



id = "california_79"
temp_temps <- temps$temperature[temps$unique.ID == id] + 3.4

preds <- calculate_fitness(temp_temps, origin = "california", return_hourly = T)
rel_diff <- ((cumsum(preds$hourly$lambda.composite) / seq_along(preds$hourly$lambda.composite))^(365*24) - (cumsum(preds$hourly$lambda.separate) / seq_along(preds$hourly$lambda.separate))^(365*24)) / ((cumsum(preds$hourly$lambda.separate) / seq_along(preds$hourly$lambda.separate))^(365*24) + (cumsum(preds$hourly$lambda.composite) / seq_along(preds$hourly$lambda.composite))^(365*24))

daily_matrix <- matrix(temp_temps, nrow = 24)
daily_means <- colMeans(daily_matrix, na.rm = TRUE)
daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)

plot(NA, 
     xlim = c(1, 365), 
     xlab = "Day of Year", 
     ylab = "Temp (°C)", 
     main = "Temperature Regime",
     ylim = c(0, 50))

polygon(x = c(1:365, 365:1),
        y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "40"),
        border = NA)
lines(1:365, daily_means[1:365], col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)      

plot(preds$hourly$lambda.composite^1~preds$hourly$hour, type = "l", xlab = "", col = "white", ylim = c(-0.1, 2.2), yaxt = "n", xaxt = "n", ylab = "Lambda")
axis(2, c(0, 1, 1.008^96), c(0, 1, 1.008))
axis(4, c(-1, 0, 1))
axis(1, c(1, 365*24), c("Jan 1", "Dec 31"))

polygon(y = c(0, preds$hourly$lambda.composite^96, 0),
        x = c(0:(length(preds$hourly$lambda.composite)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "30"))
polygon(y = c(0, preds$hourly$lambda.separate^96, 0),
        x = c(0:(length(preds$hourly$lambda.separate)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "50"))
lines(rel_diff+1, col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)


id = "california_85"
temp_temps <- temps$temperature[temps$unique.ID == id] + 4.5

preds <- calculate_fitness(temp_temps, origin = "california", return_hourly = T)
rel_diff <- (exp(cumsum(log(preds$hourly$lambda.composite)) / seq_along(preds$hourly$lambda.composite))^(365*24) - exp(cumsum(log(preds$hourly$lambda.separate)) / seq_along(preds$hourly$lambda.separate))^(365*24)) / (exp(cumsum(log(preds$hourly$lambda.composite)) / seq_along(preds$hourly$lambda.composite))^(365*24) + exp(cumsum(log(preds$hourly$lambda.separate)) / seq_along(preds$hourly$lambda.separate))^(365*24))

daily_matrix <- matrix(temp_temps, nrow = 24)
daily_means <- colMeans(daily_matrix, na.rm = TRUE)
daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)

plot(NA, 
     xlim = c(1, 365), 
     ylim = range(c(daily_mins, daily_maxs), na.rm = TRUE), 
     xlab = "Day of Year", 
     ylab = "Temp (°C)", 
     main = "Temperature Regime")

polygon(x = c(1:365, 365:1),
        y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "40"),
        border = NA)
lines(1:365, daily_means[1:365], col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)      

plot(preds$hourly$lambda.composite^1~preds$hourly$hour, type = "l", xlab = "", col = "white", ylim = c(-0.1, 2.2), yaxt = "n", xaxt = "n", ylab = "Lambda")
axis(2, c(0, 1, 1.008^96), c(0, 1, 1.008))
axis(4, c(-1, 0, 1))
axis(1, c(1, 365*24), c("Jan 1", "Dec 31"))

polygon(y = c(0, preds$hourly$lambda.composite^96, 0),
        x = c(0:(length(preds$hourly$lambda.composite)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "30"))
polygon(y = c(0, preds$hourly$lambda.separate^96, 0),
        x = c(0:(length(preds$hourly$lambda.separate)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "50"))
lines(rel_diff + 1, col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)


id = "brazil_74"
temp_temps <- temps$temperature[temps$unique.ID == id] + 1.5

preds <- calculate_fitness(temp_temps, origin = "california", return_hourly = T)
rel_diff <- ((cumsum(preds$hourly$lambda.composite) / seq_along(preds$hourly$lambda.composite))^(365*24) - (cumsum(preds$hourly$lambda.separate) / seq_along(preds$hourly$lambda.separate))^(365*24)) / ((cumsum(preds$hourly$lambda.separate) / seq_along(preds$hourly$lambda.separate))^(365*24) + (cumsum(preds$hourly$lambda.composite) / seq_along(preds$hourly$lambda.composite))^(365*24))

daily_matrix <- matrix(temp_temps, nrow = 24)
daily_means <- colMeans(daily_matrix, na.rm = TRUE)
daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)

plot(NA, 
     xlim = c(1, 365), 
     xlab = "Day of Year", 
     ylab = "Temp (°C)", 
     main = "Temperature Regime",
     ylim = c(0, 50))

polygon(x = c(1:365, 365:1),
        y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "40"),
        border = NA)
lines(1:365, daily_means[1:365], col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)      

plot(preds$hourly$lambda.composite^1~preds$hourly$hour, type = "l", xlab = "", col = "white", ylim = c(-0.1, 2.2), yaxt = "n", xaxt = "n", ylab = "Lambda")
axis(2, c(0, 1, 1.008^96), c(0, 1, 1.008))
axis(4, c(-1, 0, 1))
axis(1, c(1, 365*24), c("Jan 1", "Dec 31"))

polygon(y = c(0, preds$hourly$lambda.composite^96, 0),
        x = c(0:(length(preds$hourly$lambda.composite)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "30"))
polygon(y = c(0, preds$hourly$lambda.separate^96, 0),
        x = c(0:(length(preds$hourly$lambda.separate)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "50"))
lines(rel_diff+1, col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)



id = "brazil_22"
temp_temps <- temps$temperature[temps$unique.ID == id] + 1

preds <- calculate_fitness(temp_temps, origin = "california", return_hourly = T)
rel_diff <- (exp(cumsum(log(preds$hourly$lambda.composite)) / seq_along(preds$hourly$lambda.composite))^(365*24) - exp(cumsum(log(preds$hourly$lambda.separate)) / seq_along(preds$hourly$lambda.separate))^(365*24)) / (exp(cumsum(log(preds$hourly$lambda.composite)) / seq_along(preds$hourly$lambda.composite))^(365*24) + exp(cumsum(log(preds$hourly$lambda.separate)) / seq_along(preds$hourly$lambda.separate))^(365*24))

daily_matrix <- matrix(temp_temps, nrow = 24)
daily_means <- colMeans(daily_matrix, na.rm = TRUE)
daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)

plot(NA, 
     xlim = c(1, 365), 
     ylim = range(c(daily_mins, daily_maxs), na.rm = TRUE), 
     xlab = "Day of Year", 
     ylab = "Temp (°C)", 
     main = "Temperature Regime")

polygon(x = c(1:365, 365:1),
        y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "40"),
        border = NA)
lines(1:365, daily_means[1:365], col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)      

plot(preds$hourly$lambda.composite^1~preds$hourly$hour, type = "l", xlab = "", col = "white", ylim = c(-0.1, 2.2), yaxt = "n", xaxt = "n", ylab = "Lambda")
axis(2, c(0, 1, 1.008^96), c(0, 1, 1.008))
axis(4, c(-1, 0, 1))
axis(1, c(1, 365*24), c("Jan 1", "Dec 31"))

polygon(y = c(0, preds$hourly$lambda.composite^96, 0),
        x = c(0:(length(preds$hourly$lambda.composite)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "30"))
polygon(y = c(0, preds$hourly$lambda.separate^96, 0),
        x = c(0:(length(preds$hourly$lambda.separate)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "50"))
lines(rel_diff + 1, col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)



id = "yemen_25"
temp_temps <- temps$temperature[temps$unique.ID == id] + 1

preds <- calculate_fitness(temp_temps, origin = "california", return_hourly = T)
rel_diff <- ((cumsum(preds$hourly$lambda.composite) / seq_along(preds$hourly$lambda.composite))^(365*24) - (cumsum(preds$hourly$lambda.separate) / seq_along(preds$hourly$lambda.separate))^(365*24)) / ((cumsum(preds$hourly$lambda.separate) / seq_along(preds$hourly$lambda.separate))^(365*24) + (cumsum(preds$hourly$lambda.composite) / seq_along(preds$hourly$lambda.composite))^(365*24))

daily_matrix <- matrix(temp_temps, nrow = 24)
daily_means <- colMeans(daily_matrix, na.rm = TRUE)
daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)

plot(NA, 
     xlim = c(1, 365), 
     xlab = "Day of Year", 
     ylab = "Temp (°C)", 
     main = "Temperature Regime",
     ylim = c(0, 50))

polygon(x = c(1:365, 365:1),
        y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "40"),
        border = NA)
lines(1:365, daily_means[1:365], col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)      

plot(preds$hourly$lambda.composite^1~preds$hourly$hour, type = "l", xlab = "", col = "white", ylim = c(-0.1, 2.2), yaxt = "n", xaxt = "n", ylab = "Lambda")
axis(2, c(0, 1, 1.008^96), c(0, 1, 1.008))
axis(4, c(-1, 0, 1))
axis(1, c(1, 365*24), c("Jan 1", "Dec 31"))

polygon(y = c(0, preds$hourly$lambda.composite^96, 0),
        x = c(0:(length(preds$hourly$lambda.composite)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "30"))
polygon(y = c(0, preds$hourly$lambda.separate^96, 0),
        x = c(0:(length(preds$hourly$lambda.separate)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "50"))
lines(rel_diff+1, col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)


id = "yemen_67"
temp_temps <- temps$temperature[temps$unique.ID == id] + 1.5

preds <- calculate_fitness(temp_temps, origin = "california", return_hourly = T)
rel_diff <- (exp(cumsum(log(preds$hourly$lambda.composite)) / seq_along(preds$hourly$lambda.composite))^(365*24) - exp(cumsum(log(preds$hourly$lambda.separate)) / seq_along(preds$hourly$lambda.separate))^(365*24)) / (exp(cumsum(log(preds$hourly$lambda.composite)) / seq_along(preds$hourly$lambda.composite))^(365*24) + exp(cumsum(log(preds$hourly$lambda.separate)) / seq_along(preds$hourly$lambda.separate))^(365*24))

daily_matrix <- matrix(temp_temps, nrow = 24)
daily_means <- colMeans(daily_matrix, na.rm = TRUE)
daily_mins <- apply(daily_matrix, 2, min, na.rm = TRUE)
daily_maxs <- apply(daily_matrix, 2, max, na.rm = TRUE)

plot(NA, 
     xlim = c(1, 365), 
     ylim = range(c(daily_mins, daily_maxs), na.rm = TRUE), 
     xlab = "Day of Year", 
     ylab = "Temp (°C)", 
     main = "Temperature Regime")

polygon(x = c(1:365, 365:1),
        y = c(daily_mins[1:365], rev(daily_maxs[1:365])), 
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "40"),
        border = NA)
lines(1:365, daily_means[1:365], col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)      

plot(preds$hourly$lambda.composite^1~preds$hourly$hour, type = "l", xlab = "", col = "white", ylim = c(-0.1, 2.2), yaxt = "n", xaxt = "n", ylab = "Lambda")
axis(2, c(0, 1, 1.008^96), c(0, 1, 1.008))
axis(4, c(-1, 0, 1))
axis(1, c(1, 365*24), c("Jan 1", "Dec 31"))

polygon(y = c(0, preds$hourly$lambda.composite^96, 0),
        x = c(0:(length(preds$hourly$lambda.composite)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "30"))
polygon(y = c(0, preds$hourly$lambda.separate^96, 0),
        x = c(0:(length(preds$hourly$lambda.separate)+1)),
        border = F,
        col = paste0(unique(dat_wide$color[dat_wide$unique.ID == id]), "50"))
lines(rel_diff + 1, col = unique(dat_wide$color[dat_wide$unique.ID == id]), lwd = 2)


##########
dev.off()

pdf("./output/Fig_3_regional_fitness_diff_map_ari.pdf", width = 4, height = 3)
#### Plot map of differences between composite and separate arithmetic fitness  ----

for(warm in c(0, 2.5, 5)){
  
  # Subset data for only warming == 0
  dat_map <- dat_wide[dat_wide$warming == warm, ]
  
  # Apply the asymmetric color mapping
  point_colors <- colorRampPalette(c("darkcyan", "white", "red"))(101)[dat_map$fitness.diff.ari*50+51]
  
  xlim_ca <- c(-125, -114)
  ylim_ca <- c(32, 42)
  
  xlim_br <- c(-75, -34)
  ylim_br <- c(-35.5, 6)
  
  xlim_ye <- c(42, 55)
  ylim_ye <- c(9, 21)
  
  # Plot 1: California
  map("state", # Use the "state" database
      fill = TRUE, col = "gray80", bg = "gray95",
      xlim = xlim_ca, ylim = ylim_ca)
  title("California")
  # Add the points
  points(dat_map$longitude, dat_map$latitude,
         pch = 21,
         bg = point_colors,
         lwd = 0.5,
         cex = 1) 
  
  # Plot 2: Brazil
  map("world", 
      fill = TRUE, col = "gray80", bg = "gray95",
      xlim = xlim_br, ylim = ylim_br)
  title("Brazil")
  # Add the points
  points(dat_map$longitude, dat_map$latitude,
         pch = 21,
         bg = point_colors,
         lwd = 0.5,
         cex = 1)
  
  # Add a legend to the middle plot
  legend("bottom",
         legend = c(paste("Underest. (-1)", sep = ""),
                    "No diff. (0)",
                    paste("Overest. (1)", sep = "")),
         fill = c("darkcyan", "white", "red"),
         horiz = TRUE, 
         bty = "n",   
         cex = 1)   
  
  # Plot 3: Yemen 
  map("world", 
      fill = TRUE, col = "gray80", bg = "gray95",
      xlim = xlim_ye, ylim = ylim_ye)
  title("Yemen")
  # Add the points
  points(dat_map$longitude, dat_map$latitude,
         pch = 21,
         bg = point_colors,
         lwd = 0.5,
         cex = 1)
  
}

## Plot gradient

# Get the corresponding colors
legend_colors <- colorRampPalette(c("darkcyan", "white", "red"))(101)

# Create a horizontal matrix (1 row) of these colors
color_matrix <- as.raster(matrix(legend_colors, nrow = 1))

# Set up an empty plot. 
plot(NA, xlim = c(-1, 1), ylim = c(0, 2), 
     xaxt = 'n', yaxt = 'n', bty = 'n', xlab = "", ylab = "")

# Plot the color bar
rasterImage(color_matrix, 
            xleft = -1, 
            xright = 1, 
            ybottom = 0.5, 
            ytop = 1)


# Add numeric labels below the bar
text(x = c(-1, 0, 1), 
     y = 0.4,      
     pos = 1,        
     labels = c(-1, 0, 1),
     cex = 0.8)     

text(x = 0, y = 1.5, pos = 3, 
     labels = "Fitness Difference\n(Underestimate < 0 > Overestimate)",
     cex = 0.8,
     font = 2) 

##########
dev.off()

pdf("./output/Fig_3_regional_fitness_diff_map_geo.pdf", width = 4, height = 3)
#### Plot map of differences between composite and separate geometric fitness  ----

for(warm in c(0, 2.5, 5)){
  
  # Subset data for only warming == 0
  dat_map <- dat_wide[dat_wide$warming == warm, ]
  
  # Apply the asymmetric color mapping
  point_colors <- colorRampPalette(c("darkcyan", "white", "red"))(101)[dat_map$fitness.diff.geo*50+51]
  
  xlim_ca <- c(-125, -114)
  ylim_ca <- c(32, 42)
  
  xlim_br <- c(-75, -34)
  ylim_br <- c(-35.5, 6)
  
  xlim_ye <- c(42, 55)
  ylim_ye <- c(9, 21)
  
  # Plot 1: California
  map("state", # Use the "state" database
      fill = TRUE, col = "gray80", bg = "gray95",
      xlim = xlim_ca, ylim = ylim_ca)
  title("California")
  # Add the points
  points(dat_map$longitude, dat_map$latitude,
         pch = 21,
         bg = point_colors,
         lwd = 0.5,
         cex = 1) 
  
  # Plot 2: Brazil
  map("world", 
      fill = TRUE, col = "gray80", bg = "gray95",
      xlim = xlim_br, ylim = ylim_br)
  title("Brazil")
  # Add the points
  points(dat_map$longitude, dat_map$latitude,
         pch = 21,
         bg = point_colors,
         lwd = 0.5,
         cex = 1)
  
  # Add a legend to the middle plot
  legend("bottom",
         legend = c(paste("Underest. (-1)", sep = ""),
                    "No diff. (0)",
                    paste("Overest. (1)", sep = "")),
         fill = c("darkcyan", "white", "red"),
         horiz = TRUE, 
         bty = "n",   
         cex = 1)   
  
  # Plot 3: Yemen 
  map("world", 
      fill = TRUE, col = "gray80", bg = "gray95",
      xlim = xlim_ye, ylim = ylim_ye)
  title("Yemen")
  # Add the points
  points(dat_map$longitude, dat_map$latitude,
         pch = 21,
         bg = point_colors,
         lwd = 0.5,
         cex = 1)
  
}

## Plot gradient

# Get the corresponding colors
legend_colors <- colorRampPalette(c("darkcyan", "white", "red"))(101)

# Create a horizontal matrix (1 row) of these colors
color_matrix <- as.raster(matrix(legend_colors, nrow = 1))

# Set up an empty plot. 
plot(NA, xlim = c(-1, 1), ylim = c(0, 2), 
     xaxt = 'n', yaxt = 'n', bty = 'n', xlab = "", ylab = "")

# Plot the color bar
rasterImage(color_matrix, 
            xleft = -1, 
            xright = 1, 
            ybottom = 0.5, 
            ytop = 1)


# Add numeric labels below the bar
text(x = c(-1, 0, 1), 
     y = 0.4,      
     pos = 1,        
     labels = c(-1, 0, 1),
     cex = 0.8)     

text(x = 0, y = 1.5, pos = 3, 
     labels = "Fitness Difference\n(Underestimate < 0 > Overestimate)",
     cex = 0.8,
     font = 2) 

##########
dev.off()

pdf("./output/Fig_3_regional_fitness_diff_density_ari.pdf", width = 2.5, height = 2.5)
#### Plot densities of differences in arithmetic fitness for different warming levels / locations ----

col_palette <- colorRampPalette(c("darkcyan", "white", "red"))(101)

for(warm in c(0, 2.5, 5)){
  
  for(loc in c("california", "brazil", "yemen")){
    
    # Subset data
    dat_dens <- dat_wide[dat_wide$warming == warm & dat_wide$location == loc, ]
    
    # Handle cases with insufficient data
    if(nrow(dat_dens) < 2) {
      plot(1, type = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "",
           main = paste(loc, "at", warm, "°C\n(No Data)", sep = " "))
      next 
    }
    
    rmse_val <- sqrt(mean(dat_dens$fitness.diff.ari^2, na.rm = TRUE)) # Added ^2 for RMSE definition logic check (Root Mean SQUARE Error)
    rmse_text <- paste("RMSE:", sprintf("%.3f", rmse_val))
    
    # Calculate density
    d <- density(dat_dens$fitness.diff.ari, bw = 0.03, n = 512)
    d$y <- d$y/max(d$y)
    
    # Set up the empty plot
    plot(d, type = 'n',
         main = paste(loc, "at", warm, "°C", sep = " "),
         xlab = "", xaxt = "n",
         ylab = "", yaxt = "n",
         ylim = c(0, 1),
         xlim = c(-1, 1)) 
    
    # Calculate midpoints of the density segments
    rect_x_mid <- (d$x[-length(d$x)] + d$x[-1]) / 2
    
    # --- COLOR MAPPING FIX ---
    # Clamp values to [-1, 1] range
    vals_clamped <- pmax(-1, pmin(1, rect_x_mid))
    
    # Calculate index: -1 -> 1, 0 -> 51, 1 -> 101
    col_idx <- round(vals_clamped * 50 + 51)
    rect_colors <- col_palette[col_idx]
    # -------------------------
    
    # Draw gradient rectangles
    rect(xleft = d$x[-length(d$x)], 
         ybottom = 0,
         xright = d$x[-1],             
         ytop = d$y[-length(d$y)],    
         col = rect_colors,
         border = NA) 
    
    # Redraw density line
    lines(d, col = 'black', lwd = 1.5)
    
    legend("topleft", legend = rmse_text, bty = "n")
    abline(v = 0, lty = 2, col = 'grey30')
  }
}

##########
dev.off()

pdf("./output/Fig_3_regional_fitness_diff_density_geo.pdf", width = 2.5, height = 2.5)
#### Plot densities of differences in geometric fitness for different warming levels / locations ----

# Define the palette consistent with the map
col_palette <- colorRampPalette(c("darkcyan", "white", "red"))(101)
for(warm in c(0, 2.5, 5)){
  
  for(loc in c("california", "brazil", "yemen")){
    
    # Subset data
    dat_dens <- dat_wide[dat_wide$warming == warm & dat_wide$location == loc, ]
    
    # Handle cases with insufficient data
    if(nrow(dat_dens) < 2) {
      plot(1, type = "n", xaxt = 'n', yaxt = 'n', xlab = "", ylab = "",
           main = paste(loc, "at", warm, "°C\n(No Data)", sep = " "))
      next 
    }
    
    rmse_val <- sqrt(mean(dat_dens$fitness.diff.geo^2, na.rm = TRUE)) # Adjusted to proper RMSE calculation
    rmse_text <- paste("RMSE:", sprintf("%.3f", rmse_val))
    
    # Calculate density
    d <- density(dat_dens$fitness.diff.geo, bw = 0.03, n = 512)
    d$y <- d$y/max(d$y)
    
    # Set up empty plot
    plot(d, type = 'n',
         main = paste(loc, "at", warm, "°C", sep = " "),
         xlab = "", xaxt = "n",
         ylab = "", yaxt = "n",
         ylim = c(0, 1),
         xlim = c(-1, 1)) 
    
    # Calculate midpoints
    rect_x_mid <- (d$x[-length(d$x)] + d$x[-1]) / 2
    
    # --- COLOR MAPPING FIX ---
    # Clamp values to [-1, 1] range
    vals_clamped <- pmax(-1, pmin(1, rect_x_mid))
    
    # Calculate index: -1 -> 1, 0 -> 51, 1 -> 101
    col_idx <- round(vals_clamped * 50 + 51)
    rect_colors <- col_palette[col_idx]
    # -------------------------
    
    # Draw gradient rectangles
    rect(xleft = d$x[-length(d$x)], 
         ybottom = 0,
         xright = d$x[-1],             
         ytop = d$y[-length(d$y)],    
         col = rect_colors,
         border = NA) 
    
    # Redraw density line
    lines(d, col = 'black', lwd = 1.5)
    
    legend("topleft", legend = rmse_text, bty = "n")
    abline(v = 0, lty = 2, col = 'grey30')
  }
}
##########
dev.off()

pdf("./output/Fig_S3_fitness_diff_histograms.pdf", width = 2.5, height = 2.5)
#### Plot histograms of arithmetic fitness differences when separate component function reaches 0 or not ----

hist(dat_wide[dat_wide$geomean.fitness.composite != 0,]$fitness.diff.ari, breaks = 10, xlim = c(-1, 1), ylim = c(0, 7000))
hist(dat_wide[dat_wide$geomean.fitness.composite == 0,]$fitness.diff.ari, breaks = 50, xlim = c(-1, 1), ylim = c(0, 7000))

##########
dev.off()

#### Growth rate difference between the models across mean temperature and temperature variability ----

# Intrinsic growth rate r = log(geometric-mean lambda) of each model per site
# x warming combination, and the yearly difference between the models where
# both predict a viable population. Combinations where the composite model
# gives lambda = 0 while the separate model is viable are counted separately
# (the reverse case does not occur); both-zero combinations are excluded
r_sep <- log(dat_wide$geomean.fitness.separate)
r_comp <- log(dat_wide$geomean.fitness.composite)
both_viable <- dat_wide$geomean.fitness.separate > 0 & dat_wide$geomean.fitness.composite > 0
both_zero <- dat_wide$geomean.fitness.separate == 0 & dat_wide$geomean.fitness.composite == 0
comp_zero_only <- dat_wide$geomean.fitness.composite == 0 & dat_wide$geomean.fitness.separate > 0
r_diff_year <- ifelse(both_viable, 365 * 24 * (r_comp - r_sep), NA)

# Mean and SD of the hourly temperatures of the year for which fitness is
# simulated, with the site mean shifted by warming
mean_temps <- aggregate(temperature ~ unique.ID, data = temps[temps$middle.year == 1, ], "mean")
sd_temps <- aggregate(temperature ~ unique.ID, data = temps[temps$middle.year == 1, ], "sd")
thermal_regimes <- merge(mean_temps, sd_temps, by = "unique.ID")
colnames(thermal_regimes) <- c("unique.ID", "mean_temp", "sd_temp")
mean_temp <- thermal_regimes$mean_temp[match(dat_wide$unique.ID, thermal_regimes$unique.ID)] + dat_wide$warming
sd_temp <- thermal_regimes$sd_temp[match(dat_wide$unique.ID, thermal_regimes$unique.ID)]

site_colors <- c(california = "#1f78b4", brazil = "#33a02c", yemen = "#ff7f00")

# Fill colors: yearly r difference (composite - separate), clipped to +/- log 2
# (a two-fold difference in yearly growth potential); navy: composite lambda = 0
# while the separate model is viable
diff_palette <- colorRampPalette(c("darkcyan", "white", "red"))(101)
diff_cap <- log(2)
diff_color <- function(r_diff) diff_palette[round(pmin(pmax(r_diff, -diff_cap), diff_cap) / diff_cap * 50 + 51)]
diff_bar_labels <- c(expression("" <= -ln ~ 2), "0", expression("" >= ln ~ 2))
diff_bar_title <- expression(italic(r)[comp] - italic(r)[sep] ~ (yr^-1))
collapse_col <- "navy"

# Color bar and grey swatches in the top margin, positioned relative to the plot frame
draw_diff_legend <- function(bar_palette, bar_labels, bar_title, swatch_col, swatch_label) {

  usr <- par("usr")
  x_span <- diff(usr[1:2])
  bar_y <- usr[4] + c(0.06, 0.11) * diff(usr[3:4])

  bar_x <- seq(usr[1] + 0.03 * x_span, usr[1] + 0.27 * x_span, length.out = 102)
  rect(head(bar_x, -1), bar_y[1], tail(bar_x, -1), bar_y[2], col = bar_palette, border = NA, xpd = NA)
  rect(bar_x[1], bar_y[1], bar_x[102], bar_y[2], border = "grey50", lwd = 0.5, xpd = NA)
  text(c(bar_x[1], mean(bar_x), bar_x[102]), bar_y[2], bar_labels, adj = c(0.5, -0.4), cex = 0.6, xpd = NA)
  text(bar_x[102] + 0.02 * x_span, mean(bar_y), bar_title, adj = c(0, 0.5), cex = 0.6, xpd = NA)

  swatch_x <- usr[2] - rev(seq_along(swatch_col)) * 0.19 * x_span
  rect(swatch_x, bar_y[1], swatch_x + 0.03 * x_span, bar_y[2], col = swatch_col, border = "grey50", lwd = 0.5, xpd = NA)
  text(swatch_x + 0.04 * x_span, mean(bar_y), swatch_label, adj = c(0, 0.5), cex = 0.6, xpd = NA)
}

##########

pdf("./output/Fig_S5_growth_rate_diff.pdf", width = 5, height = 4.5)
#### Plot every site as a point in mean x SD space at three warming levels ----

# Sites at present-day, +2.5 and +5 C in one panel, filled by the yearly r
# difference (composite - separate); navy: composite lambda = 0 while the
# separate model is viable; light grey: both lambda = 0; squares: opposite
# signs of r. Drawing order: both-zero, composite collapse, then increasing
# |difference| so the most different sites end up on top, squares last
shown <- dat_wide$warming %in% c(0, 2.5, 5)
point_category <- ifelse(both_zero, 0, ifelse(comp_zero_only, 1, 2))
sign_flip <- both_viable & (r_sep > 0) != (r_comp > 0)
point_col <- c("#e5e5e5", collapse_col, NA)[point_category + 1]
point_col[point_category == 2] <- diff_color(r_diff_year[point_category == 2])
point_rank <- ifelse(point_category == 2, abs(r_diff_year), 0)
point_order <- which(shown)[order(point_category[shown], sign_flip[shown], point_rank[shown])]

plot(mean_temp[point_order], sd_temp[point_order],
     pch = ifelse(sign_flip[point_order], 22, 21), cex = ifelse(sign_flip[point_order], 0.9, 0.7), lwd = 0.3,
     col = ifelse(sign_flip[point_order], "black", "grey70"), bg = point_col[point_order],
     xlab = "Mean (°C)", ylab = "SD (°C)")
draw_diff_legend(diff_palette, diff_bar_labels, diff_bar_title,
                 c(collapse_col, "#e5e5e5"), c(expression(lambda[comp] == 0), expression(both ~ lambda == 0)))
legend("bottomleft", legend = expression(opposite ~ signs ~ of ~ italic(r)), pch = 22, pt.bg = "white",
       col = "black", pt.cex = 0.9, bty = "n", cex = 0.7)

##########
dev.off()

#### Time to a two-fold divergence between the models ----

# Starting from the same population size, the two predictions differ by a
# factor of two after log(2) / |r_comp - r_sep| years, with the composite
# prediction ending up higher (r_comp > r_sep) or lower (r_comp < r_sep).
# Combinations where the composite model gives lambda = 0 while the separate
# model is viable diverge immediately (t = 0, composite lower); both-zero
# combinations never diverge and are excluded.
# Growth regimes: all combinations; lambda < 1 in exactly one model (which
# includes every composite collapse of a growing population); lambda < 1 in
# both models; lambda > 1 in both models

time_to_twofold <- ifelse(both_viable, log(2) / abs(r_diff_year), NA)
time_to_twofold[comp_zero_only] <- 0
diverging <- !both_zero
composite_higher <- both_viable & r_diff_year > 0
composite_lower <- comp_zero_only | (both_viable & r_diff_year < 0)

growth_regime <- list("all" = diverging,
                      "lambda < 1 in one model" = diverging & xor(r_sep < 0, r_comp < 0),
                      "lambda < 1 in both models" = diverging & r_sep < 0 & r_comp < 0,
                      "lambda > 1 in both models" = both_viable & r_sep > 0 & r_comp > 0)

time_marks <- c(1 / 365, 1 / 12, 1, 10)

# Share of the selected combinations diverged by each time, in one direction
cumulative_share <- function(sel, direction, at) {
  sapply(at, function(x) sum(sel & direction & time_to_twofold <= x)) / sum(sel) * 100
}

twofold_summary <- do.call(rbind, lapply(names(growth_regime), function(regime) {
  do.call(rbind, lapply(c("brazil", "california", "yemen", "all"), function(loc) {
    sel <- growth_regime[[regime]] & (loc == "all" | dat_wide$location == loc)
    higher <- cumulative_share(sel, composite_higher, time_marks)
    lower <- cumulative_share(sel, composite_lower, time_marks)
    data.frame(regime = regime,
               location = loc,
               n = sum(sel),
               pct.1.day = higher[1] + lower[1],
               pct.1.month = higher[2] + lower[2],
               pct.1.year = higher[3] + lower[3],
               pct.10.years = higher[4] + lower[4],
               pct.composite.higher.1.year = higher[3],
               pct.composite.lower.1.year = lower[3],
               pct.composite.higher.10.years = higher[4],
               pct.composite.lower.10.years = lower[4],
               median.years = median(time_to_twofold[sel]))
  }))
}))
print(twofold_summary)

write.table(twofold_summary, "./data/time_to_twofold_summary.txt", sep = "\t", row.names = F, quote = F)

##########

pdf("./output/Fig_3_time_to_twofold.pdf", width = 9, height = 2.5)
#### Plot the cumulative share of combinations that have diverged two-fold ----

# One panel per growth regime, one line per region and direction: composite
# prediction higher above the axis, lower below. Points mark the shares at
# 1 day, 1 month, 1 year and 10 years, labelled above the points; dotted
# lines at 0, 50 and 100 %. Curves that never leave zero, and points at
# zero, are omitted
time_grid <- 10^seq(log10(1 / 365), 1, length.out = 400)
regions <- c("brazil", "yemen", "california")
regime_titles <- list(expression(All ~ combinations),
                      expression(lambda < 1 ~ "in one model"),
                      expression(lambda < 1 ~ "in both models"),
                      expression(lambda > 1 ~ "in both models"))

# Keep coinciding labels at least 16 percentage points apart, pushing upward and staying below 100
dodge_labels <- function(y) {
  o <- order(y)
  y_sorted <- y[o]
  for (k in seq_along(y_sorted)[-1]) y_sorted[k] <- max(y_sorted[k], y_sorted[k - 1] + 16, na.rm = TRUE)
  y[o] <- y_sorted - max(0, max(y_sorted, na.rm = TRUE) - 100)
  y
}

par(mfrow = c(1, 4))
for (i in seq_along(growth_regime)) {

  plot(NA, xlim = range(time_grid), ylim = c(-100, 100), log = "x", xaxt = "n", yaxt = "n",
       xlab = "Time", ylab = "", main = regime_titles[[i]])
  axis(1, at = time_marks, labels = c("1 d", "1 mo", "1 yr", "10 yr"))
  axis(2, at = seq(-100, 100, by = 50), labels = c("100", "50", "0", "50", "100"))
  mtext(c("Composite higher", "Separate higher"), side = 2, at = c(65, -65), line = 2.5, cex = 0.5)
  abline(h = seq(-100, 100, by = 50), lty = 3, col = "gray60")

  for (direction in c("higher", "lower")) {
    direction_sign <- ifelse(direction == "higher", 1, -1)
    direction_sel <- if (direction == "higher") composite_higher else composite_lower
    share_marks <- sapply(regions, function(loc) {
      cumulative_share(growth_regime[[i]] & dat_wide$location == loc, direction_sel, time_marks)
    })
    label_y <- t(apply(direction_sign * share_marks, 1, dodge_labels))
    label_pos <- c(4, 3, 3, 3)

    for (loc in regions) {
      sel <- growth_regime[[i]] & dat_wide$location == loc
      if (!any(sel)) next
      share_grid <- cumulative_share(sel, direction_sel, time_grid)
      if (!any(share_grid > 0)) next
      lines(time_grid, direction_sign * share_grid, col = site_colors[loc], lwd = 1)
      marked <- share_marks[, loc] > 0
      points(time_marks[marked], direction_sign * share_marks[marked, loc], pch = 16, col = site_colors[loc], cex = 0.8)
      if (any(marked)) {
        text(time_marks[marked], label_y[marked, loc], paste0(round(share_marks[marked, loc]), "%"),
             pos = label_pos[marked], offset = 0.4, cex = 0.7, col = site_colors[loc], xpd = NA)
      }
    }
  }
}
par(mfrow = c(1, 1))

##########
dev.off()
