#### Preparations ----

# Read data
dat <- read.delim("./data/thermal_performance.txt")
dat$nonviable <- ifelse(dat$adult.offspring == 0, 1, 0)
dat$temperature.factor <- factor(dat$temperature)
dat$dummy <- ifelse(dat$replicate == "brazil1", 0.5, ifelse(dat$replicate == "brazil2", -0.5, 0))
dat$color <- ifelse(dat$origin == "ca", "#1f78b4", ifelse(dat$origin == "yem", "#ff7f00", "#33a02c"))

# Packages
library(brms)
library(bayestestR)

# Load functions
source("./R/functions.R")

##########

#### Model lifetime reproductive success ----

# Specify priors
LRS_priors <- prior(normal(29, 5), nlpar="ToptFec", lb = 0) +
  prior(normal(70,25), nlpar="fecmax", lb = 0) +
  prior(lognormal(log(22), log(1.1)), nlpar="Tb", lb = 0) +
  prior(lognormal((2.5),(1.5)), nlpar="K", lb = 0) + 
  prior(lognormal(log(0.03), log(1.2)), nlpar="a", lb = 0) +
  prior(normal(28, 2), nlpar="ToptViab") +
  prior(normal(-4.5, 2), nlpar="k")+
  prior(normal(0, 1), nlpar = "rep") + 
  prior(normal(0, 2.5), dpar="shape") 

# Initialize starting values
init <- parse(text=paste0("list(b_ToptFec = as.array(rnorm(3, 29, 1.5)),
  b_fecmax = as.array(rnorm(3, 70, 1.5)),
  b_Tb = as.array(rlnorm(3, log(25), log(1.1))),
  b_K = as.array(rlnorm(1, log(4), log(1.5)))
  )"))
set.seed(1646);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))


LRS_model <- brm(bf(adult.offspring ~ log(fecundity(temperature, ToptFec, fecmax, Tb, K)) + dummy * rep,
                          nlf(zi ~ nonviability(temperature, a, ToptViab, k)),
                          shape ~ 0 + temperature.factor,
                          K + a + ToptViab + k + rep ~ 1,
                          ToptFec + fecmax + Tb ~ 0 + origin,
                          nl = T),      
                       data = dat, 
                       family = zero_inflated_negbinomial,
                       prior = LRS_priors,
                       control = list(adapt_delta = 0.95, max_treedepth = 15),
                       warmup = 2000,
                       init = inits_list,
                       stanvar = stan_funs,
                       iter = 3000,
                       thin = 1,
                       seed = 5557,
                       cores = 6,
                       chains = 6,
                       file = "./models/LRS_model")

pdf("./output/Fig_S7_pp_check_LRS.pdf", width = 10, height = 3)
pp_check(LRS_model, type = "intervals_grouped", group = "origin")
dev.off()
summary(LRS_model)
#plot(LRS_model)

##########

#### Save LRS model parameter estimates ----

LRS_model_posteriors <- as.data.frame(LRS_model)

## Extract posterior estimates and uncertainties
# Initialize an empty data frame to store the results
posterior_estimates <- data.frame(trait=character(),
                                  parameter=character(),
                                  effect=character(),
                                  origin=character(),
                                  MAP=numeric(),
                                  hpdi.lower=numeric(),
                                  hpdi.upper=numeric(),
                                  stringsAsFactors=FALSE)

# Parameters in the development rate model
parameters <- c("ToptFec", "fecmax", "Tb", "K", "a", "ToptViab", "k", "rep")
trait <- "fecundity"

# Main effects: parameters depending on origin and selection_regime
for (parameter in parameters) {
  for (origin in unique(dat$origin)) {
    
    effect <- "main"
    if (parameter %in% c("K", "a", "ToptViab", "k", "rep")){
      param_name <- paste0("b_", parameter, "_Intercept")
  }
    else{
    param_name <- paste0("b_", parameter, "_origin", origin)
    }
    if (param_name %in% names(LRS_model_posteriors)) {
      samples <- LRS_model_posteriors[[param_name]]
      # Compute MAP and HPDI
      MAP <- point_estimate(samples, centrality='MAP')
      HPDI <- ci(samples, method='HDI', ci=0.9)
      
      posterior_estimates <- rbind(posterior_estimates,
                                   data.frame(trait=trait,
                                              parameter=parameter,
                                              effect=effect,
                                              origin=origin,
                                              MAP=MAP,
                                              hpdi.lower=HPDI$CI_low,
                                              hpdi.upper=HPDI$CI_high,
                                              stringsAsFactors=FALSE))
    } else {
      # Parameter not found
      message(paste("Parameter", param_name, "not found in posterior samples."))
    }
  }
}

# Rename the posterior estimates data frame
posterior_estimates_LRS <- posterior_estimates
write.table(posterior_estimates_LRS, file = "./data/parameter_estimates_LRS.txt", sep = "\t", row.names = FALSE, quote = F)

##########

pdf("./output/Fig_S2_fitted_LRS_reaction_norms.pdf", width = 2.5, height = 3)
#### Plot fecundity and non-viability ----

# Generate dataframe for predictions
temp_seq <- seq(10, 45, by = 0.25)
newdata <- expand.grid(
  temperature = temp_seq,
  origin = unique(dat$origin),
  temperature.factor = "29",
  dummy = 0
)

# Get fitted values for the fecundity component
fitted_vals <- fitted(LRS_model, newdata = newdata, re_formula = NULL, dpar = "mu")

# Combine
preds <- cbind(newdata, fitted_vals)

# Plot
plot(NA,
     ylim = c(0, 100),
     xlim = c(10, 45),
     xlab = "Temperature (°C)", ylab = "Adult offspring")
points(dat$adult.offspring ~ jitter(dat$temperature, 2), col = dat$color, pch = 16, cex = 0.5)

# Plot loop
for(origin in unique(dat$origin)){
  
  # Subset data 
  dat_rep <- dat[dat$origin == origin, ]
  
  # Subset predictions
  pred_rep <- preds[preds$origin == origin, ]
  
  # Plot raw data
  
  # 95% CI
  polygon(
    x = c(pred_rep$temperature, rev(pred_rep$temperature)),
    y = (c(pred_rep$Q2.5, rev(pred_rep$Q97.5))),
    col = paste0(unique(dat$color[dat$origin == origin]), "30"),
    border = NA
  )
  
  # Add the fitted line (the posterior mean of the curve)
  lines(Estimate ~ temperature, data = pred_rep, col = unique(dat$color[dat$origin == origin]), lwd = 2)
}


# Generate dataframe for predictions
temp_seq <- seq(10, 45, by = 0.25)
newdata <- expand.grid(
  temperature = temp_seq,
  origin = unique(dat$origin)[1],
  temperature.factor = "29",
  dummy = 0
)

# Get fitted values for the fecundity component
fitted_vals_zi <- fitted(LRS_model, newdata = newdata, dpar = "zi", re_formula = NA)

# Combine
preds_zi <- cbind(newdata, fitted_vals_zi)

# 95% CI
polygon(
  x = c(preds_zi$temperature, rev(preds_zi$temperature)),
  y = (c(preds_zi$Q2.5, rev(preds_zi$Q97.5)))*100,
  col = "#00000030",
  border = NA
)

# Add the fitted line (the posterior mean of the curve)
lines(Estimate*100 ~ temperature, data = preds_zi, col = "#000000", lwd = 3)

##########
dev.off()

#### Model development rate ----

# Specify priors
devrate_priors <- prior(normal(10, 3), nlpar = "Tmin", lb = 0) +
  prior(normal(32.5,3), nlpar = "Topt") +
  prior(normal(37.5, 3), nlpar = "Tmax") +
  prior(uniform(0, 1), nlpar = "rmax", lb = 0, ub = 1) +
  prior(normal(0, 1), nlpar = "rep") 

# Initialize starting values
init <- parse(text=paste0("list(b_Tmin = as.array(rnorm(3, 10, 1)),
  b_Topt = as.array(rnorm(3, 30, 1)),
  b_Tmax = as.array(rnorm(3, 40, 1)),
  b_rmax = as.array(rnorm(3, 0.5, 0.01))
  )"))
set.seed(1646);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))


devrate_model <- brm(bf(dev.rate ~ log(LRF(temperature, Tmin, Topt, Tmax, rmax)) + rep * dummy,
                    Tmin + Topt + Tmax + rmax ~ 0 + origin,
                    rep ~ 1,
                    nl = T),      
                 data = dat, 
                 family = lognormal,
                 prior = devrate_priors,
                 control = list(adapt_delta = 0.95, max_treedepth = 15),
                 warmup = 2000,
                 init = inits_list,
                 stanvar = stan_funs,
                 iter = 3000,
                 thin = 1,
                 seed = 5559,
                 cores = 6,
                 chains = 6,
                 file = "./models/devrate_model")

pdf("./output/Fig_S7_pp_check_dev.pdf", width = 10, height = 3)
pp_check(devrate_model, type = "intervals_grouped", group = "origin")
dev.off()
summary(devrate_model)
#plot(devrate_model)
#plot(conditional_effects(devrate_model))

##########

pdf("./output/Fig_S2_fitted_development_reaction_norms.pdf", width = 2.5, height = 3)
#### Plot development rate  ----

LRF <- Vectorize(LRF)
# Generate dataframe for predictions
temp_seq <- seq(10, 45, by = 0.25)
newdata <- expand.grid(
  temperature = temp_seq,
  origin = unique(dat$origin),
  dummy = 0
)

# Get fitted values for the fecundity component
fitted_vals <- fitted(devrate_model, newdata = newdata)

# Combine
preds <- cbind(newdata, fitted_vals)

plot(NA,
     ylim = c(0, 1 / (24 * 15)),
     xlim = c(10, 45),
     xlab = "Temperature (°C)", ylab = "Adult offspring")
points(dat$dev.rate ~ jitter(dat$temperature, 2), col = dat$color, pch = 16, cex = 0.5)


# Plot loop
for(origin in unique(dat$origin)){
  
  # Subset data 
  dat_rep <- dat[dat$origin == origin, ]
  
  # Subset predictions
  pred_rep <- preds[preds$origin == origin, ]
  
  # 95 % CI
  polygon(
    x = c(pred_rep$temperature, rev(pred_rep$temperature)),
    y = (c(pred_rep$Q2.5, rev(pred_rep$Q97.5))),
    col = paste0(unique(dat$color[dat$origin == origin]), "30"),
    border = NA
  )
  
  # Add the fitted line (the posterior mean of the curve)
  lines(Estimate ~ temperature, data = pred_rep, col = unique(dat$color[dat$origin == origin]), lwd = 3)
}


##########
dev.off()

#### Save dev. rate posterior estimates ----
dev_rate_model_posteriors <- as.data.frame(devrate_model)

## Extract posterior estimates and uncertainties
# Initialize an empty data frame to store the results
posterior_estimates <- data.frame(trait=character(),
                                  parameter=character(),
                                  effect=character(),
                                  origin=character(),
                                  MAP=numeric(),
                                  hpdi.lower=numeric(),
                                  hpdi.upper=numeric(),
                                  stringsAsFactors=FALSE)

# Parameters in the development rate model
parameters <- c("Tmin", "Topt", "Tmax", "rmax", "rep")
trait <- "development_rate"

# Main effects: parameters depending on origin and selection_regime
for (parameter in parameters) {
  for (origin in unique(dat$origin)) {

      effect <- "main"
      param_name <- paste0("b_", parameter, "_origin", origin)
      if (parameter == "rep"){
        param_name <- paste0("b_", parameter, "_Intercept")
      }
      else{
        param_name <- paste0("b_", parameter, "_origin", origin)
      }
      
      if (param_name %in% names(dev_rate_model_posteriors)) {
        samples <- dev_rate_model_posteriors[[param_name]]
        # Compute MAP and HPDI
        MAP <- point_estimate(samples, centrality='MAP')
        HPDI <- ci(samples, method='HDI', ci=0.9)
        
        posterior_estimates <- rbind(posterior_estimates,
                                     data.frame(trait=trait,
                                                parameter=parameter,
                                                effect=effect,
                                                origin=origin,
                                                MAP=MAP,
                                                hpdi.lower=HPDI$CI_low,
                                                hpdi.upper=HPDI$CI_high,
                                                stringsAsFactors=FALSE))
      } else {
        # Parameter not found
        message(paste("Parameter", param_name, "not found in posterior samples."))
      }
  }
}

# Rename the posterior estimates data frame
posterior_estimates_dev_rate <- posterior_estimates
write.table(posterior_estimates_dev_rate, file = "./data/parameter_estimates_dev_rate.txt", sep = "\t", row.names = FALSE, quote = F)

##########

#### Model growth rate ----

parameter_estimates_dev_rate <- subset(read.delim("./data/parameter_estimates_dev_rate.txt"), effect == "main" & parameter == "Tmax" | effect == "main" & parameter == "Tmin" | effect == "main" & parameter == "Topt")
specific_priors_list <- list()
for (i in 1:nrow(parameter_estimates_dev_rate)) {
  
  nlpar_name <- parameter_estimates_dev_rate$parameter[i]
  origin_val <- parameter_estimates_dev_rate$origin[i]
  map_estimate <- parameter_estimates_dev_rate$MAP[i]
  prior_sd <- 1 
  coef_name <- paste0("b_", nlpar_name, "_origin", origin_val)
  
  current_prior <- eval(parse(text = paste0("prior(normal(", map_estimate, ", ", prior_sd, "),
                         nlpar = ", nlpar_name, ",
                         coef = ", gsub("b_.*?_", "", coef_name), ")")))
  
  specific_priors_list[[length(specific_priors_list) + 1]] <- current_prior
}

# Combine the generated priors into a single brmsprior object
specific_priors_list <- Reduce(`+`, specific_priors_list)

# Specify priors
growthrate_priors <-  prior(uniform(0, 1), nlpar = "rmax", lb = 0, ub = 1) +
  prior(normal(0, 1), nlpar = "rep") +
  specific_priors_list

# Initialize starting values
init <- parse(text=paste0("list(b_Tmin = as.array(rnorm(3, 10, 1)),
  b_Topt = as.array(rnorm(3, 30, 1)),
  b_Tmax = as.array(rnorm(3, 40, 1)),
  b_rmax = as.array(rnorm(3, 0.5, 0.01))
  )"))

set.seed(6093);inits_list <- list(eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init), eval(init))

growthrate_model <- brm(bf(growth.rate ~ log(LRF(temperature, Tmin, Topt, Tmax, rmax)) + dummy * rep,
                        Tmin + Topt + Tmax + rmax ~ 0 + origin,
                        rep ~ 1,
                        nl = T),      
                     data = dat, 
                     family = lognormal,
                     prior = growthrate_priors,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     warmup = 2000,
                     init = inits_list,
                     stanvar = stan_funs,
                     iter = 3000,
                     thin = 1,
                     seed = 1893,
                     cores = 6,
                     chains = 6,
                     file = "./models/growthrate_model")
pdf("./output/Fig_S7_pp_check_growth.pdf", width = 10, height = 3)
pp_check(growthrate_model, type = "intervals_grouped", group = "origin")
dev.off()
summary(growthrate_model)
#plot(growthrate_model)

##########

pdf("./output/Fig_S2_fitted_growth_reaction_norms.pdf", width = 2.5, height = 3)
#### Plot growth rate  ----

LRF <- Vectorize(LRF)
# Generate dataframe for predictions
temp_seq <- seq(10, 45, by = 0.25)
newdata <- expand.grid(
  temperature = temp_seq,
  origin = unique(dat$origin),
  dummy = 0
)

# Get fitted values for the fecundity component
fitted_vals <- fitted(growthrate_model, newdata = newdata)

# Combine
preds <- cbind(newdata, fitted_vals)

plot(NA,
     ylim = c(0, 0.000005),
     xlim = c(10, 45),
     xlab = "Temperature (°C)", ylab = "Growth rate (g/h)")
points(dat$growth.rate ~ jitter(dat$temperature, 2), col = dat$color, pch = 16, cex = 0.5)

# Plot loop
for(origin in unique(dat$origin)){
  
  # Subset data 
  dat_rep <- dat[dat$origin == origin, ]
  
  # Subset predictions
  pred_rep <- preds[preds$origin == origin, ]
  
  # 95% CI
  polygon(
    x = c(pred_rep$temperature, rev(pred_rep$temperature)),
    y = (c(pred_rep$Q2.5, rev(pred_rep$Q97.5))),
    col = paste0(unique(dat$color[dat$origin == origin]), "30"),
    border = NA
  )
  
  # Add the fitted line (the posterior mean of the curve)
  lines(Estimate ~ temperature, data = pred_rep, col = unique(dat$color[dat$origin == origin]), lwd = 3)
}


##########
dev.off()

#### Save growth rate posterior estimates ----

growth_rate_model_posteriors <- as.data.frame(growthrate_model)

## Extract posterior estimates and uncertainties
# Initialize an empty data frame to store the results
posterior_estimates <- data.frame(trait=character(),
                                  parameter=character(),
                                  effect=character(),
                                  origin=character(),
                                  MAP=numeric(),
                                  hpdi.lower=numeric(),
                                  hpdi.upper=numeric(),
                                  stringsAsFactors=FALSE)

# Parameters in the development rate model
parameters <- c("Tmin", "Topt", "Tmax", "rmax", "rep")
trait <- "growth_rate"

# Main effects: parameters depending on origin and selection_regime
for (parameter in parameters) {
  for (origin in unique(dat$origin)) {
    
    effect <- "main"
    param_name <- paste0("b_", parameter, "_origin", origin)
    if (parameter == "rep"){
      param_name <- paste0("b_", parameter, "_Intercept")
    }
    else{
      param_name <- paste0("b_", parameter, "_origin", origin)
    }
    
    if (param_name %in% names(growth_rate_model_posteriors)) {
      samples <- growth_rate_model_posteriors[[param_name]]
      # Compute MAP and HPDI
      MAP <- point_estimate(samples, centrality='MAP')
      HPDI <- ci(samples, method='HDI', ci=0.9)
      
      posterior_estimates <- rbind(posterior_estimates,
                                   data.frame(trait=trait,
                                              parameter=parameter,
                                              effect=effect,
                                              origin=origin,
                                              MAP=MAP,
                                              hpdi.lower=HPDI$CI_low,
                                              hpdi.upper=HPDI$CI_high,
                                              stringsAsFactors=FALSE))
    } else {
      # Parameter not found
      message(paste("Parameter", param_name, "not found in posterior samples."))
    }
  }
}

# Rename the posterior estimates data frame
posterior_estimates_growth_rate <- posterior_estimates
write.table(posterior_estimates_growth_rate, file = "./data/parameter_estimates_growth_rate.txt", sep = "\t", row.names = FALSE, quote = F)

##########
