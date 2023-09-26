# pdf plotting to manually check SDM outputs for irregularities

library(raster)
library(sf)
library(rnaturalearth)

# creating land
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
land_poly <- as_Spatial(land)

# plotting
pdf('./plots/all_sdms.pdf')
for(i in seq_along(list.files('./final_sdm_outputs'))) {
  want_ras <- raster(paste0('./final_sdm_outputs/', list.files('./final_sdm_outputs')[i]))
  spec_name <- gsub('_', ' ', names(want_ras))
  plot(land_poly, col = NA, main = paste0("SDM for ", spec_name))
  plot(want_ras, add = T)
}
dev.off()