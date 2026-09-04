#### Preparations ----

# Read data
temps <- read.delim("./data/era5_thermal_regimes.txt")

library(maps)

# Load functions
source("./R/functions.R")

##########

pdf("./output/Fig_2_map.pdf", width = 5, height = 3.5)
#### Plot map ----


unique_coords <- unique(temps[, c("longitude", "latitude", "location")])
countries <- unique(unique_coords$location)
country_colors <- setNames(c("#33a02c", "#1f78b4", "#ff7f00"), countries)
point_colors <- country_colors[unique_coords$location]
xlims <- range(unique_coords$longitude) + c(-5, 5)
ylims <- range(unique_coords$latitude) + c(-5, 5)


map("world", fill = TRUE, col = "gray80", bg = "gray95", 
    mar = c(0, 0, 2, 0), xlim = xlims, ylim = ylims)

# Add the points to the map
points(unique_coords$longitude, unique_coords$latitude, 
       pch = 21,
       bg = point_colors, 
       lwd = 0.5,
       cex = 0.7)         

##########
dev.off()

length(unique(temps$unique.ID))

