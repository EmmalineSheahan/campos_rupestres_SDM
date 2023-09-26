# making raster stacks and richness maps

library(sf)
library(terra)
library(sp)
library(rnaturalearth)
library(geobr)
library(rgdal)
library(raster)
library(dplyr)

# creating new CRS
moll_crs <- "+proj=moll +lon_0=1 +datum=WGS84 +units=m +no_defs"

# creating land
land <- ne_countries(country = "Brazil", returnclass = "sf", scale = 50)
land_poly <- as_Spatial(land)
land_poly <- spTransform(land_poly, CRSobj = moll_crs)

# creating phytogeographic regions
biomes <- read_biomes(year = 2019, simplified = TRUE)
biomes_no_sea <- biomes[-7,]

biomes_st <- st_as_sf((biomes_no_sea), remove = FALSE, crs = 4326, agr = "constant")[1]
biomes_st$name_biome <- c("Amazonia", "Caatinga", "Cerrado", "Mata_Atlantica", "Pampa",
                          "Pantanal")

cr <- shapefile('./rasters/cr.shp')

first_shp <- raster::union(as_Spatial(biomes_st %>% filter(name_biome == "Caatinga")), 
                   as_Spatial(biomes_st %>% filter(name_biome == "Cerrado")))
regions <- raster::union(first_shp, 
                         as_Spatial(biomes_st %>% filter(name_biome == "Mata_Atlantica")))

plot(cr, border = "red")
plot(regions, col = NA, add = T)

regions_mol <- spTransform(regions, CRSobj = moll_crs)
writeOGR(regions_mol, dsn = './rasters/', layer = 'regions_mol',
         driver = 'ESRI Shapefile', overwrite = T)

cerrado <- as_Spatial(biomes_st %>% filter(name_biome == "Cerrado"))
cerrado <- spTransform(cerrado, CRSobj = moll_crs)

caatinga <- as_Spatial(biomes_st %>% filter(name_biome == "Caatinga"))
caatinga <- spTransform(caatinga, CRSobj = moll_crs)

mata <- as_Spatial(biomes_st %>% filter(name_biome == "Mata_Atlantica"))
mata <- spTransform(mata, CRSobj = moll_crs)

# creating a base raster to transform each sdm to
base_raster <- raster(ext = extent(regions_mol), crs = moll_crs, resolution = 1000)

# species list
# move all of the tifs under hulls to the final_sdm_output directory
all_species <- list.files('./final_sdm_outputs')
all_species <- gsub('.tif', '', all_species)

# creating directory for transformed rasters
dir.create('./transformed_rasters')

# transforming and cropping rasters to region of interest, writing to transformed_rasters dir
for (i in seq_along(all_species)) {
  test_raster <- raster(paste0('./final_sdm_outputs/', all_species[i], '.tif'))
  test_points <- rasterToPoints(test_raster)
  test_points <- data.frame(test_points)
  colnames(test_points) <- c("x", "y", "z")
  test_points_sp <- test_points
  coordinates(test_points_sp) <- ~x+y
  proj4string(test_points_sp) <- "+proj=longlat +ellps=WGS84 +no_defs"
  test_points_moll <- spTransform(test_points_sp, CRSobj = moll_crs)
  new_raster <- rasterize(test_points_moll, base_raster, 
                            field = test_points$z)
  new_raster <- mask(new_raster, regions_mol)
  names(new_raster) <- all_species[i]
  writeRaster(new_raster, filename = paste0('./transformed_rasters/', 
                                              all_species[i], '.tif'),
                format = "GTiff", overwrite = T)
}

# create list of filenames for stack
species_list <- list.files('./transformed_rasters')
for (i in 1:length(species_list)) {
  species_list[i] <- paste0('./transformed_rasters/', species_list[i])
}
species_list <- as.list(species_list)

# create stack
species_stack <- stack(species_list)
species_spatrast <- rast(species_stack)
writeRaster(species_spatrast, filename = './data/all_species_stack.tif',
            overwrite = T)

# create richness map
cr_richness <- calc(species_stack, fun = sum, na.rm = T)
cr_richness <- mask(cr_richness, regions_mol)
writeRaster(cr_richness, filename = './data/campos_rupestres_richness.tif',
            format = "GTiff", overwrite = T)

pdf('./plots/richness_campos_rupestres.pdf')
plot(land_poly, col = NA, main = paste0('Species Richness in the Campos Rupestres'))
plot(cr_richness, add = T)
dev.off()

pdf('./plots/richness_campos_rupestres_phytoregions.pdf')
plot(cr_richness, main = paste0('Species Richness in the Campos Rupestres'))
plot(land_poly, col = NA, add = T)
plot(cerrado, border = "blue", col = NA, add = T)
plot(caatinga, border = "red", col = NA, add = T)
plot(mata, border = "purple", col = NA, add = T)
dev.off()
