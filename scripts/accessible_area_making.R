# creating accessible area polygons per species

library(sp)
library(rangeBuilder)
library(alphahull)
library(sf)
library(dplyr)
library(rgeos)
library(geodata)
library(ENMeval)

# creating species list
species_list <- list.files('./Occurrences')
species_names <- gsub('.csv', '', species_list)

# checking how many records exist for each species
records_amount <- vector(length = length(species_list))
 for (i in seq_along(species_list)) {
   t <- read.csv(paste0('./Occurrences/', species_list[i]))
   records_amount[i] <- nrow(t)
 }

# we'll have to drop species with too few records
drop_species <- species_names[which(records_amount < 4)]
amount_of_records_drop <- records_amount[which(records_amount < 4)]
drop_species_df <- data.frame(drop_species, amount_of_records_drop)
colnames(drop_species_df) <- c("Species", "Amount_of_Records")
write.csv(drop_species_df, file = './data/species_with_under_4_coords.csv')

# for species with 4-7 points, we can make alpha hulls around them for richness and PD
# these will have to be vetted by hand
species_list_new <- species_list[-which(records_amount < 4)]
species_names_new <- species_names[-which(records_amount < 4)]
records_amount_new <- vector(length = length(species_list_new))
for (i in seq_along(species_list_new)) {
   t <- read.csv(paste0('./Occurrences/', species_list_new[i]))
   records_amount_new[i] <- nrow(t)
}
 
hull_species <- species_names_new[which(records_amount_new < 8)]
records_hull <- records_amount_new[which(records_amount_new < 8)]
hull_species_df <- data.frame(hull_species, records_hull)
colnames(hull_species_df) <- c("Species", "Amount_of_Records")
write.csv(hull_species_df, file = './data/hull_species_with_4_to_7_coords.csv')

# writing the final list of species with 8 or greater cleaned occurrences for modeling
species_list_final <- species_list_new[-which(records_amount_new < 8)]
species_names_final <- species_names_new[-which(records_amount_new < 8)]
species_titles_final <- gsub('_', ' ', species_names_final)
write.table(species_names_final, file = './data/modeled_species_list.txt', row.names = F,
             col.names = F)

# creating land
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
land_poly <- as_Spatial(land)
land_poly_cea <- spTransform(land_poly, CRSobj = "+proj=cea +lat_ts=0 +lon_0")

#creating directory for accessible area shapefiles
dir.create('./Acc_Area')

# creating various accessible areas of buffer distances 2, 3, 4, and 5 times larger
# than the area of the coordinates, and then running initial models to select the best
# buffer for each species via AIC
buff_props <- c(2, 3, 4, 5)

# pulling down environmental rasters from WorldClim
dir.create('./rasters')
all_clim <- worldclim_global(var= "bio", res = 5, path = './rasters')
names(all_clim) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
                     "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", 
                     "bio16", "bio17", "bio18", "bio19")
wanted_clim <- c(all_clim[[1]], all_clim[[4]], all_clim[[7]], all_clim[[12]], 
                 all_clim[[13]], all_clim[[14]], all_clim[[15]])

# creating accessible area polygons
# for some reason there's an inexplicable error occurring in getDynamicAlphaHull
# even though there are enough coordinates for the function to work,
# i'm inserting a tryCatch for now which will store the species it's refusing
# to work on

alphahull_malfunction <- integer()

pdf('./plots/accessible_areas_modeled_poales.pdf')
for (i in seq_along(species_list_final)) {
  tryCatch(
    {
  # creating the different buffers
  temp_spec_mat <- read.csv(paste0('./Occurrences/', species_list_final[i]))
  temp_spec_mat <- temp_spec_mat[,2:4]
  temp_spec_coords <- temp_spec_mat[,1:2]
  temp_spec <- temp_spec_mat
  coordinates(temp_spec) <- ~Longitude+Latitude
  proj4string(temp_spec) <- CRS("+proj=longlat +datum=WGS84")
  temp_trans <- spTransform(temp_spec, CRSobj = "+proj=cea +lat_ts=0 +lon_0")
  temp_alph <- getDynamicAlphaHull(temp_spec_coords, fraction = 1, partCount = 1, 
                                   initialAlpha = 20, clipToCoast = "terrestrial",
                                   coordHeaders = c("Longitude", "Latitude"))
  temp_alph <- spTransform(as_Spatial(temp_alph[[1]]), 
                           CRSobj = "+proj=cea +lat_ts=0 +lon_0")
  shape_list <- vector("list", length = length(buff_props))
  for (j in seq_along(buff_props)) {
    buffDist = (sqrt(buff_props[j]*gArea(temp_alph)) - sqrt(gArea(temp_alph)))/2
    shape_new <- raster::buffer(x = temp_alph, width = buffDist, dissolve = T)
    shape_new_clipped <- intersect(shape_new, land_poly_cea)
    if(gArea(shape_new_clipped) > (gArea(land_poly_cea)*0.8)) {
      shape_list[[j]] <- NA
    } else {
    shape_sf <- st_as_sf(shape_new_clipped)
    shape_final <- st_transform(shape_sf, "+proj=longlat +datum=WGS84")
    if(!(any(st_is_valid(shape_final)))) {
      shape_final <- st_make_valid(shape_final)
    }
    shape_list[[j]] <- shape_final
    }
  }
  
  # j loop to crop rasters to each polygon, run ENMeval for each polygon, and store
  # model results for each polygon
  model_results <- vector(length = length(shape_list))
  for (j in 1:length(shape_list)) {
    if(is.logical(shape_list[[j]])) {
      model_results[j] <- NA
    } else {
    cropped_clim <- crop(wanted_clim, shape_list[[j]])
    clipped_clim <- mask(cropped_clim, shape_list[[j]])
    simple_model <- ENMevaluate(occs = temp_spec_mat[,1:2], envs = clipped_clim, 
                                algorithm = "maxnet",
                                tune.args = list(fc = c("L","LQ"), rm = 1:2),
                                partitions = "none", n.bg = 10000)
    mod_res <- min(eval.results(simple_model)$AICc)
    model_results[j] <- mod_res
    }
  }
  
  # find which polygon produced the lowest AIC
  want_shape_num <- which(model_results == min(model_results, na.rm = T))
  want_shape <- shape_list[[want_shape_num]]
  
  # create plots to ensure accessible areas are accurate
  plot(land_poly, col = NA, main = paste0("Accessible Area for ", 
                                          species_titles_final[i]))
  plot(as_Spatial(want_shape), col = "yellow", add = T)
  plot(temp_spec, col = "red", pch = 19, cex = 0.25, add = T)
  
  # write that polygon to the ./Acc_Area directory
  st_write(want_shape, dsn = './Acc_Area', layer = species_names_final[i], 
           driver = "ESRI Shapefile",
           append = F)
    }, 
    error=function(err){
    message('On iteration ',i, ' there was an error: ',err)
    alphahull_malfunction <<-c(alphahull_malfunction,i)
  }
  )
}
dev.off()

malfunctioned_species <- species_list_final[alphahull_malfunction]
write.table(malfunctioned_species, file = './data/malfunctioned_species_list.txt', 
            row.names = F,
            col.names = F)

