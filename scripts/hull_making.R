# creating alpha hulls for species with 3-7 points

library(dplyr)
library(sf)
library(terra)
library(rangeBuilder)
library(sp)
library(raster)

# reading in lists
hull_species1 <- read.csv('./data/hull_species_with_4_to_7_coords.csv')
hull_species2 <- read.csv('./data/species_with_under_4_coords.csv')

hull_species2 <- hull_species2 %>% filter(Amount_of_Records == 3)
wanted_species <- c(hull_species1$Species, hull_species2$Species)
wanted_species_title <- gsub('_', ' ', wanted_species)

# creating land
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
land_poly <- as_Spatial(land)

# import campos rupestres polygon
camposrupestres <- shapefile('./rasters/cr.shp')
proj4string(camposrupestres) <- CRS("+proj=longlat +datum=WGS84")

# creating hull directory
dir.create('./hulls')

# loop to create final occurrence rasters for hull species
pdf('./plots/hull_species.pdf')
hull_malfunction <- integer()
for(i in seq_along(wanted_species)) {
  tryCatch(
    {
  occs <- read.csv(paste0('./Occurrences/', wanted_species[i], '.csv'))
  coordinates(occs) <- ~Longitude+Latitude
  proj4string(occs) <- CRS("+proj=longlat +datum=WGS84")
  hull_poly <- getDynamicAlphaHull(x = occs@coords, fraction = 1, clipToCoast = "terrestrial")
  whole_rast <- rast(crs = "+proj=longlat +datum=WGS84", resolution = 0.008333333, 
                     extent = ext(st_bbox(hull_poly[[1]])), vals = 1)
  clipped_rast <- mask(whole_rast, st_as_sf(hull_poly[[1]]))
  plot(camposrupestres, col = NA, border = "red", main = wanted_species_title[i])
  plot(land_poly, col = NA, add = T)
  plot(raster(clipped_rast), add = T)
  names(clipped_rast) <- wanted_species[i]
  writeRaster(clipped_rast, filename = paste0('./hulls/', 
                                          wanted_species[i], '.tif'), filetype = "GTiff",
              overwrite = T)
    }, 
  error=function(err){
    message('On iteration ',task_id, ' there was an error: ',err)
    hull_malfunction <<-c(hull_malfunction,task_id)
  }
  )
}
dev.off()
