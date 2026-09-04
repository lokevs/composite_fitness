#### Preparations ----

# Read data
dat <- read.delim("./data/thermal_performance.txt")
dat$nonviable <- ifelse(dat$adult.offspring == 0, 1, 0)
dat$temperature.factor <- factor(dat$temperature)

# Load functions
source("./R/functions.R")

##########

pdf("./output/Fig_S1_conceptual_figure_components.pdf", width = 2.5, height = 3)
#### Plots of composite / ind. components  ----

# Define a bunch of made-up reaction norms
devrate_slow <- function(x) LRF(x, 10, 33, 38, 1 / (24 * 2))
growthrate_slow <- function(x) LRF(x, 8, 32, 40, 1 / (24 * 2))
fec_slow <- function(x, mass) fecundity(x, 26, 1 * (mass / 1.5), 23, 2)
viab_slow <- function(x) fert_rate(x, 10, 33, 38, 1 / (24 * 2), 0.025, 26, -2.5)
fitness_composite_slow <- function(x) {
  devtime <- 1 / devrate_slow(x)
  mass <- growthrate_slow(x) * devtime
  offspring <- (viab_slow(x) ^ devtime) * fec_slow(x, mass)
  return(ifelse(is.nan(offspring / devtime), 0, offspring / devtime))
}
fitness_components_slow <- function(x){
  if(max(devrate_slow(x)) <= 0){
    return(0)
  }
  else{
    temps <- x[1:which(cumsum(devrate_slow(x)) > 1)[1]]
    mass <- sum(growthrate_slow(temps))
    offspring <- mean(fec_slow(temps, mass))
    fitness <- offspring * prod(viab_slow(temps)) / length(temps)
    return(fitness)
  }
}

# Define a bunch of made-up reaction norms
devrate_fast <- function(x) LRF(x, 7, 33, 42, 0.03125000)
growthrate_fast <- function(x) LRF(x, 7, 32, 43, 0.03279963)
fec_fast <- function(x, mass) fecundity(x, 25, 0.7 * (mass / 1.8), 35, 3)
viab_fast <- function(x) fert_rate(x, 7, 33, 42, 0.03125000, 0.05221412 , 26, -2.5)
fitness_composite_fast <- function(x) {
  devtime <- 1 / devrate_fast(x)
  mass <- growthrate_fast(x) * devtime
  offspring <- (viab_fast(x) ^ devtime) * fec_fast(x, mass)
  return(ifelse(is.nan(offspring / devtime), 0, offspring / devtime))
}
fitness_components_fast <- function(x){
  if(max(devrate_fast(x)) <= 0){
    return(0)
  }
  else{
    temps <- x[1:which(cumsum(devrate_fast(x)) > 1)[1]]
    mass <- sum(growthrate_fast(temps))
    offspring <- mean(fec_fast(temps, mass))
    fitness <- offspring * prod(viab_fast(temps)) / length(temps)
    return(fitness)
  }
}


# Plot dev rate
plot(y = devrate_slow(5:45) / max(devrate_fast(5:45)), x = 5:45,
     ylim = c(0, 1.2),
     type = "l",
      lwd = 2,
      ylab = "", yaxt = "n",
      xlab = "Temperature (°C)")
lines(y = devrate_fast(5:45) / max(devrate_fast(5:45)), x = 5:45,
      lwd = 2,
      lty = 1, 
      col = "orangered")

# Plot growth rate
plot(y = growthrate_slow(5:45) / max(growthrate_fast(5:45)), x = 5:45,
     ylim = c(0, 1.2),
     type = "l",
     lwd = 2,
     ylab = "", yaxt = "n",
     xlab = "Temperature (°C)")
lines(y = growthrate_fast(5:45) / max(growthrate_fast(5:45)), x = 5:45,
      lwd = 2,
      lty = 1, 
      col = "orangered")

# Plot viability rate
plot(y = viab_slow(5:45), x = 5:45,
     ylim = c(0, 1),
     type = "l",
     lwd = 2,
     ylab = "", yaxt = "n",
     xlab = "Temperature (°C)")
lines(y = viab_fast(5:45), x = 5:45,
      lwd = 2,
      lty = 1, 
      col = "orangered")

# Plot fecundity
plot(y = fec_slow(5:45, mass = 1) / max(fec_slow(5:45, mass = 1)), x = 5:45,
     ylim = c(0, 1.2),
     type = "l",
     lwd = 2,
     ylab = "", yaxt = "n",
     xlab = "Temperature (°C)")
lines(y = fec_fast(5:45, mass = 1)  / max(fec_slow(5:45, mass = 1)), x = 5:45,
      lwd = 2,
      lty = 1, 
      col = "orangered")


plot(y = fitness_composite_slow(5:45) / max(fitness_composite_slow(5:45), na.rm = T), x = 5:45,
     ylim = c(0, 1.2),
     type = "l",
     lwd = 2,
     ylab = "", yaxt = "n",
     xlab = "Temperature (°C)")

lines(y = fitness_composite_fast(5:45) / max(fitness_composite_slow(5:45), na.rm = T), x = 5:45,
      lwd = 2,
      lty = 1, 
      col = "orangered")


legend("topleft", title = "Life history",
       c("Slow", "Fast"),
       lty = c(1, 2), 
       col = c("black", "orangered"), 
       lwd = 2,
       cex = 1,
       bty = "n")

##########
dev.off()

pdf("./output/Fig_S1_conceptual_figure_predictions.pdf", width = 2.5, height = 3)
#### Plot predictions under fluctuations ----

temps <- 10:37
results_df <- data.frame()

# This loop can be slow, as it simulates many time steps
for(amp in c(0, 4, 8, 12, 16)){
  for(temp_mean in temps){
    # Create a long vector of fluctuating temperatures
    temporary_temps <- temp_mean + amp * sin(2 * pi * (1:(24 * 1000)) / 24)
    
    # Slow life history
    fitness_comp_slow <- fitness_components_slow(temporary_temps)
    fitness_composite_slow_val <- mean(fitness_composite_slow(temporary_temps))
    
    # Fast life history
    fitness_comp_fast <- fitness_components_fast(temporary_temps)
    fitness_composite_fast_val <- mean(fitness_composite_fast(temporary_temps))
    
    results_df <- rbind(results_df, data.frame(
      mean.temp = temp_mean,
      amplitude = amp,
      strategy = "slow",
      fitness.components = fitness_comp_slow,
      fitness.composite = fitness_composite_slow_val
    ))
    
    results_df <- rbind(results_df, data.frame(
      mean.temp = temp_mean,
      amplitude = amp,
      strategy = "fast",
      fitness.components = fitness_comp_fast,
      fitness.composite = fitness_composite_fast_val
    ))
  }
}

palettes <- list(
  grey = colorRampPalette(c("grey60", "lightgrey")),
  orange = colorRampPalette(c("#ffb399", "#ffeadb"))
)

for(amp in c(0, 4, 8, 12, 16)){

# Set up an empty plot
plot(NA,
     ylim = c(0, 1),
     xlim = c(5, 45),
     ylab = "",
     yaxt = "n",
     xlab = "Temperature (°C)")

# --- Gradient Polygon for Slow Strategy, Amplitude 10 ---
poly_data_slow <- subset(results_df, strategy == "slow" & amplitude == amp)

# Create functions that describe the top and bottom edges of the polygon
top_edge_fun <- approxfun(poly_data_slow$mean.temp, poly_data_slow$fitness.composite, rule = 2)
bottom_edge_fun <- approxfun(poly_data_slow$mean.temp, poly_data_slow$fitness.components, rule = 2)

# Define gradient properties
n_vertical_slices <- 100 # Number of vertical trapezoids for smoothness
n_gradient_steps <- 20  # Number of stacked colors per trapezoid
gradient_colors <- palettes$grey(n_gradient_steps) # Using the grey palette

x_interp <- seq(min(poly_data_slow$mean.temp), max(poly_data_slow$mean.temp), length.out = n_vertical_slices)
y_top_interp <- top_edge_fun(x_interp)
y_bottom_interp <- bottom_edge_fun(x_interp)

# Loop through each vertical slice to build the gradient
for (i in 1:(n_vertical_slices - 1)) {
  # Get the four corners of the main trapezoid for this slice
  x1 <- x_interp[i]
  x2 <- x_interp[i+1]
  y_bottom1 <- y_bottom_interp[i]
  y_bottom2 <- y_bottom_interp[i+1]
  y_top1 <- y_top_interp[i]
  y_top2 <- y_top_interp[i+1]
  
  # Skip if there's no height or data
  if (anyNA(c(y_bottom1, y_bottom2, y_top1, y_top2))) next
  if ((y_top1 - y_bottom1) <= 1e-9 && (y_top2 - y_bottom2) <= 1e-9) next
  
  # Loop vertically *within* this trapezoid to create the gradient effect
  for (j in 0:(n_gradient_steps - 1)) {
    # p1 and p2 are the proportions of the total height for this small colored polygon
    p1 <- j / n_gradient_steps      
    p2 <- (j + 1) / n_gradient_steps
    
    # Interpolate the y-coordinates for the small polygon based on the proportions
    poly_y1 <- y_bottom1 + (y_top1 - y_bottom1) * p1
    poly_y2 <- y_bottom2 + (y_top2 - y_bottom2) * p1
    poly_y3 <- y_bottom2 + (y_top2 - y_bottom2) * p2
    poly_y4 <- y_bottom1 + (y_top1 - y_bottom1) * p2
    
    #polygon(
    #  x = c(x1, x2, x2, x1),
    #  y = c(poly_y1, poly_y2, poly_y3, poly_y4),
    #  col = gradient_colors[j + 1],
    #  border = NA
    #)
  }
}

# --- Add Outlines and Other Data ---
# Add lines for the slow strategy
lines(fitness.composite/max(fitness.composite) ~ mean.temp, data = subset(results_df, strategy == "slow" & amplitude == amp)[order(subset(results_df, strategy == "slow" & amplitude == amp)$mean.temp),], lwd = 2, lty = 2, col = "grey40")
lines(fitness.components/max(fitness.composite) ~ mean.temp, data = subset(results_df, strategy == "slow" & amplitude == amp)[order(subset(results_df, strategy == "slow" & amplitude == amp)$mean.temp),], lwd = 2, lty = 1, col = "black")
}


for(amp in c(0, 4, 8, 12, 16)){
  
  # Set up an empty plot
  plot(NA,
       ylim = c(0, 1),
       xlim = c(5, 45),
       ylab = "",
       yaxt = "n",
       xlab = "Temperature (°C)")
  
  # --- Gradient Polygon for Slow Strategy, Amplitude 10 ---
  poly_data_slow <- subset(results_df, strategy == "fast" & amplitude == amp)
  
  # Create functions that describe the top and bottom edges of the polygon
  top_edge_fun <- approxfun(poly_data_slow$mean.temp, poly_data_slow$fitness.composite, rule = 2)
  bottom_edge_fun <- approxfun(poly_data_slow$mean.temp, poly_data_slow$fitness.components, rule = 2)
  
  # Define gradient properties
  n_vertical_slices <- 100 # Number of vertical trapezoids for smoothness
  n_gradient_steps <- 20  # Number of stacked colors per trapezoid
  gradient_colors <- palettes$orange(n_gradient_steps) # Using the grey palette
  
  x_interp <- seq(min(poly_data_slow$mean.temp), max(poly_data_slow$mean.temp), length.out = n_vertical_slices)
  y_top_interp <- top_edge_fun(x_interp)
  y_bottom_interp <- bottom_edge_fun(x_interp)
  
  # Loop through each vertical slice to build the gradient
  for (i in 1:(n_vertical_slices - 1)) {
    # Get the four corners of the main trapezoid for this slice
    x1 <- x_interp[i]
    x2 <- x_interp[i+1]
    y_bottom1 <- y_bottom_interp[i]
    y_bottom2 <- y_bottom_interp[i+1]
    y_top1 <- y_top_interp[i]
    y_top2 <- y_top_interp[i+1]
    
    # Skip if there's no height or data
    if (anyNA(c(y_bottom1, y_bottom2, y_top1, y_top2))) next
    if ((y_top1 - y_bottom1) <= 1e-9 && (y_top2 - y_bottom2) <= 1e-9) next
    
    # Loop vertically *within* this trapezoid to create the gradient effect
    for (j in 0:(n_gradient_steps - 1)) {
      # p1 and p2 are the proportions of the total height for this small colored polygon
      p1 <- j / n_gradient_steps      
      p2 <- (j + 1) / n_gradient_steps
      
      # Interpolate the y-coordinates for the small polygon based on the proportions
      poly_y1 <- y_bottom1 + (y_top1 - y_bottom1) * p1
      poly_y2 <- y_bottom2 + (y_top2 - y_bottom2) * p1
      poly_y3 <- y_bottom2 + (y_top2 - y_bottom2) * p2
      poly_y4 <- y_bottom1 + (y_top1 - y_bottom1) * p2
      
      #polygon(
      #  x = c(x1, x2, x2, x1),
      #  y = c(poly_y1, poly_y2, poly_y3, poly_y4),
      #  col = gradient_colors[j + 1],
      #  border = NA
      #)
    }
  }
  
  # --- Add Outlines and Other Data ---
  # Add lines for the slow strategy
  lines(fitness.composite/max(fitness.composite) ~ mean.temp, data = subset(results_df, strategy == "fast" & amplitude == amp)[order(subset(results_df, strategy == "fast" & amplitude == amp)$mean.temp),], lwd = 2, lty = 2, col = "#ffb399")
  lines(fitness.components/max(fitness.composite) ~ mean.temp, data = subset(results_df, strategy == "fast" & amplitude == amp)[order(subset(results_df, strategy == "fast" & amplitude == amp)$mean.temp),], lwd = 2, lty = 1, col = "orangered")
}


##########
dev.off()

pdf("./output/Fig_1_intro_figure_1.pdf", width = 3, height = 3.5)
#### Plot 2-process fitness (two additive) part 1 ----

# Process 1 (accumulation resource 1)
proc1 <- Vectorize(function(x){
  return(max(x-5, 0))
})

# Process 2 (accumulation resource 2)
proc2 <- Vectorize(function(x){
  return(max(10-x, 0))
})
# Plot reaction norms
curve(proc1(x), 2.5, 12.5, lty = 1, col = "grey10", lwd = 2)
curve(proc2(x), 2.5, 12.5, add = T, col = "orange", lty = 1, lwd = 2)
curve(proc2(x) * proc1(x), 2.5, 12.5, add = T, col = "red", lty = 2, lwd = 2)

E1 = 4.5
E2 = 10.5

# Variable environment
abline(v = c(E1, E2), lwd = 1, col = "grey", lty = 2)

##########
dev.off()

pdf("./output/Fig_1_intro_figure_2.pdf", width = 4, height = 3)
#### Plot 2-process fitness (two additive) part 2 ----

# Plot resource accumulation
plot(y = c(E1, E1, E2, E2, E1, E1, E2, E2),
     x = c(0, 1, 1, 2, 2, 3, 3, 4), 
     type = "l",
     lwd = 1,
     col = "grey",
     lty = 2,
     ylim = c(0, 12), 
     xlim = c(0, 6),
     ylab = "", yaxt = "n",
     xlab = "", xaxt = "n", bty = "l")
axis(2, c(4.5, 10.5), labels = c("E1", "E2"))
axis(1, c(0.5, 1.5, 2.5, 3.5), labels = c("t1", "t2", "t3", "t4"))

# Accumulation of resource 1
env <- c(E1, E2, E1, E2)
R1_vec <- c(0, cumsum(proc1(env)))
lines(0:4, R1_vec - 0.1, col = "grey10", type = "l", lwd = 2)
points(y = tail(R1_vec, 1) - 0.1, x = 4,  bg = "grey10", pch = 21, cex = 1.5)

# Accumulation of resource 2
R2_vec <- c(0, cumsum(proc2(env)))
lines(0:4, R2_vec + 0.1, col = "orange", type = "l", lwd = 2)
points(y = tail(R2_vec, 1) + 0.1, x = 4,  bg = "orange", pch = 21, cex = 1.5)

# Accumulation of composite fitness prediction (Pink dashed line)
comp_vec <- c(0, cumsum(proc1(env) * proc2(env)))
lines(0:4, comp_vec, col = "red", type = "l", lwd = 2, lty = 2)
points(y = tail(comp_vec, 1), x = 4, bg = "red", pch = 21, cex = 1.5)

# True fitness (Separate)
W_true <- tail(R1_vec, 1) * tail(R2_vec, 1)
rect(4.6, -1, 5.4, W_true / 13, col = "dodgerblue", border = NA)

# Predicted fitness (Composite)
W_pred <- sum(proc1(env) * proc2(env))
rect(5.6, -1, 6.4, W_pred / 13, col = "red", border = NA)

# Dotted baseline
curve(x * 0, 4.5, 6.5, add = TRUE, lty = 2)

##########
dev.off()

pdf("./output/Fig_1_intro_figure_3.pdf", width = 3, height = 3.5)
#### Plot 2-process fitness (one additive one multiplicative) part 1 ----

# Create file for plot

# Process 1 (mortality)
proc1 <- Vectorize(function(x){
  return(1 - exp(-5 + 0.37 * x))
})

# Process 2 (growth)
proc2 <- function(x){
  LRF(x, 0, 10, 12.5, 1)
}
curve(proc1(x), 2.5, 12.5, lty = 1, col = "grey10", lwd = 2, ylim = c(0, 1))
curve(proc2(x), 2.5, 12.5, add = T, col = "orange", lty = 1, lwd = 2)
composite <- function(x) {
  
  (proc1(x))^4 * 4*proc2(x)
  
} # 4 days life span

curve(composite(x) / 4, 2.5, 12.5, add = T, col = "red", lty = 2, lwd = 2)

E1 = 8
E2 = 12

abline(v = c(E1, E2), col = "grey", lwd = 1, lty = 2)


########
dev.off()

pdf("./output/Fig_1_intro_figure_4.pdf", width = 4, height = 3)
#### Plot 2-process fitness (one additive one multiplicative) part 2 ----

env <- c(E1, E2, E1, E2)

# Plotting parameters
plot(y = c(4.5, 4.5, 10.5, 10.5, 4.5, 4.5, 10.5, 10.5),
     x = c(0, 1, 1, 2, 2, 3, 3, 4), type = "l",
     lwd = 1,
     col = "grey",
     lty = 2,
     ylim = c(0, 12), xlim = c(0, 6),
     ylab = "", yaxt = "n", xlab = "", xaxt = "n", bty = "l")
axis(2, c(4.5, 10.5), labels = c("E1", "E2"))
axis(1, c(0.5, 1.5, 2.5, 3.5), labels = c("t1", "t2", "t3", "t4"))

# Accumulation 1: Survival (Multiplicative)
# Starts at env probability (proc1[1]), not 1
S_vec <- c(1, cumprod(proc1(env)))
lines(0:4, S_vec * 12, col = "grey10", type = "l", lwd = 2) 
points(y = tail(S_vec, 1) * 12, x = 4,  bg = "grey10", pch = 21, cex = 1.5)

# Accumulation 2: Growth (Additive)
# Starts at 0
G_vec <- c(0, cumsum(proc2(env)))
lines(0:4, G_vec * 4.6, col = "orange", type = "l", lwd = 2)
points(y = tail(G_vec, 1) * 4.6, x = 4,  bg = "orange", pch = 21, cex = 1.5)

# Accumulation of composite fitness prediction (Pink dashed line)
comp_vec <- c(0, cumsum(composite(env)) / length(env))
lines(0:4, comp_vec * 8, col = "red", type = "l", lwd = 2, lty = 2)
points(y = tail(comp_vec, 1) * 8, x = 4, bg = "red", pch = 21, cex = 1.5)

# True fitness (Final Survival * Final Mass)
W_true <- tail(S_vec, 1) * tail(G_vec, 1)
rect(4.6, -1, 5.4, W_true * 8, col = "dodgerblue", border = NA)

# Predicted fitness (Composite)
W_pred <- mean(composite(env))
rect(5.6, -1, 6.4, W_pred * 8, col = "red", border = NA)

# Dotted baseline
curve(x * 0, 4.5, 6.5, add = TRUE, lty = 2)

##########
dev.off()

pdf("./output/Fig_1_intro_figure_5.pdf", width = 3, height = 3.5)
#### Plot 2-process fitness (two additive but temporally separated) part 1 ----

# Process 1 (development t0 to t2)
proc1 <- Vectorize(function(x){
  LRF(12-x, 0, 7, 9, 1)
})

# Process 2 (growth)
proc2 <- function(x){
  LRF(x, 3, 10, 12, 1)
}
curve(proc1(x) / 1.4, 2.5, 12.5, lty = 1, col = "orange", lwd = 2, ylim = c(0, 1))
curve(proc2(x) / 1.4, 2.5, 12.5, add = T, col = "grey10", lty = 1, lwd = 2)
curve((proc2(x) + proc1(x)) / 1.4, 2.5, 12.5, add = T, col = "red", lty = 2, lwd = 2)

E1 = 5
E2 = 10
abline(v = c(E1, E2), col = "grey", lwd = 1, lty = 2)

##########
dev.off()

pdf("./output/Fig_1_intro_figure_6.pdf", width = 4, height = 3)
#### Plot 2-process fitness (two additive but temporally separated) part 2 ----

# Plotting parameters (Good conditions)
plot(y = c(4.5, 4.5, 4.5, 4.5, 10.5, 10.5, 10.5, 10.5),
     x = c(0, 1, 1, 2, 2, 3, 3, 4), type = "l",
     lwd = 1,
     col = "grey",
     lty = 2,
     ylim = c(0, 12), xlim = c(0, 6),
     ylab = "", yaxt = "n", xlab = "", xaxt = "n", bty = "l")
axis(2, c(4.5, 10.5), labels = c("E1", "E2"))
axis(1, c(0.5, 1.5, 2.5, 3.5), labels = c("t1", "t2", "t3", "t4"))

# Accumulation of development ages t0 to t2
env <- c(E1, E1, E2, E2)
R1_vec <- c(0, cumsum(proc1(env[1:2])) * 2.5)
lines(0:2, R1_vec, col = "grey10", type = "l", lwd = 2)
points(y = tail(R1_vec, 1), x = 2,  bg = "grey10", pch = 21, cex = 1.5)

# Accumulation of development ages t2 to t4
R2_vec <- c(tail(R1_vec, 1), cumsum(proc2(env[3:4])) * 2.5 + tail(R1_vec, 1)) 
lines(2:4, R2_vec , col = "orange", type = "l", lwd = 2)
points(y = tail(R2_vec, 1), x = 4,  bg = "orange", pch = 21, cex = 1.5)

# Accumulation of composite fitness prediction
comp_vec <- c(0, cumsum(proc1(env) + proc2(env)) / 2)
lines(0:4, comp_vec, col = "red", type = "l", lwd = 2, lty = 2)
points(y = tail(comp_vec, 1), x = 4, bg = "red", pch = 21, cex = 1.5)

# True fitness (Separate)
W_true <- tail(R2_vec, 1)
rect(4.6, -1, 5.4, W_true, col = "dodgerblue", border = NA)

# Predicted fitness (Composite)
W_pred <- sum(proc1(env) + proc2(env)) / 2
rect(5.6, -1, 6.4, W_pred, col = "red", border = NA)

# Dotted baseline
curve(x * 0, 4.5, 6.5, add = TRUE, lty = 2)


E1 = 10
E2 = 5

# Plotting parameters (Bad conditions)
plot(y = c(10.5, 10.5, 10.5, 10.5, 4.5, 4.5, 4.5, 4.5),
     x = c(0, 1, 1, 2, 2, 3, 3, 4), type = "l",
     lwd = 1,
     col = "grey",
     lty = 2,
     ylim = c(0, 12), xlim = c(0, 6),
     ylab = "", yaxt = "n", xlab = "", xaxt = "n", bty = "l")
axis(2, c(4.5, 10.5), labels = c("E1", "E2"))
axis(1, c(0.5, 1.5, 2.5, 3.5), labels = c("t1", "t2", "t3", "t4"))

# Accumulation of development ages t0 to t2
env <- c(E1, E1, E2, E2)
R1_vec <- c(0, cumsum(proc1(env[1:2])) * 2.5)
lines(0:2, R1_vec, col = "grey10", type = "l", lwd = 2)
points(y = tail(R1_vec, 1), x = 2,  bg = "grey10", pch = 21, cex = 1.5)

# Accumulation of development ages t2 to t4
R2_vec <- c(tail(R1_vec, 1), cumsum(proc2(env[3:4])) * 2.5 + tail(R1_vec, 1)) 
lines(2:4, R2_vec , col = "orange", type = "l", lwd = 2)
points(y = tail(R2_vec, 1), x = 4,  bg = "orange", pch = 21, cex = 1.5)

# Accumulation of composite fitness prediction
comp_vec <- c(0, cumsum(proc1(env) + proc2(env)) / 2)
lines(0:4, comp_vec, col = "red", type = "l", lwd = 2, lty = 2)
points(y = tail(comp_vec, 1), x = 4, bg = "red", pch = 21, cex = 1.5)

# True fitness (Separate)
W_true <- tail(R2_vec, 1)
rect(4.6, -1, 5.4, W_true, col = "dodgerblue", border = NA)

# Predicted fitness (Composite)
W_pred <- sum(proc1(env) + proc2(env)) / 2
rect(5.6, -1, 6.4, W_pred, col = "red", border = NA)

# Dotted baseline
curve(x * 0, 4.5, 6.5, add = TRUE, lty = 2)

##########
dev.off()
