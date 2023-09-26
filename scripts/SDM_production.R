# primary script for running ENMs and outputting final projections

library(sf)
library(terra)
library(dplyr)
library(dismo)
library(ENMeval)
library(geodata)
library(usdm)
library(blockCV)

# creating land for plotting
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
land_poly <- as_Spatial(land)
target_crs <- "+proj=longlat +datum=WGS84"

# creating directory for final outputs
dir.create('./final_sdm_outputs')
dir.create('./error_species')
dir.create('./sdm_parameters')

# environmental rasters
names_list <- c("bio1", "bio4", "bio7", "bio12", "bio13", "bio14", "bio15")
all_clim <- vector("list", length = length(names_list))
file_list <- list.files('./rasters/final_rasters')
for (i in 1:length(file_list)) {
  clim1 <- rast(paste0('./rasters/final_rasters/', file_list[i]))
  names(clim1) <- names_list[i]
  all_clim[[i]] <- clim1
}

wanted_clim <- c(all_clim[[1]], all_clim[[2]], all_clim[[3]], all_clim[[4]], 
                 all_clim[[5]], all_clim[[6]], all_clim[[7]])

# getting species lists
model_species <- read.table('./data/modeled_species_list.txt')
model_species <- model_species$V1
model_titles <- gsub('_', ' ', model_species)

acc_list <- list.files('./Acc_Area')
acc_list <- gsub('.shp', '', acc_list)
acc_list <- gsub('.shx', '', acc_list)
acc_list <- gsub('.prj', '', acc_list)
acc_list <- gsub('.dbf', '', acc_list)
acc_list <- unique(acc_list)

# trying to catch any inconsistency between the acc area files and the modeled species list
inconsistent <- which(!(model_species %in% acc_list))
if (length(inconsistent) != 0) {
  no_acc_area <- model_species[inconsistent]
  write.table(no_acc_area, file = './data/species_which_had_no_acc_area_for_modeling.txt', 
            row.names = F,
            col.names = F)
  model_species <- model_species[-inconsistent]
}

# var_remove function for single iterative variable removal
# wanted_envs = the clipped environmental spatraster used in the model
# wanted_occs = the two column occurrence matrix
# wanted_bg = the two column background matrix
# returns the spatraster with the least important variable discarded
var_remove <- function(wanted_envs, wanted_occs, wanted_bg) {
  rast_list <- vector("list", length = dim(wanted_envs)[3])
  for (j in 1:dim(wanted_envs)[3]) {
    env_new <- raster(wanted_envs[[j]])
    rast_list[[j]] <- env_new
  }
  initial_envs <- stack(rast_list)
  initial_mod <- maxent(initial_envs, p = wanted_occs, a = wanted_bg)
  initial_results <- data.frame(initial_mod@results)
  measure_name <- row.names(initial_results)
  perm_import <- data.frame(measure_name, 
                          initial_results[,1])
  colnames(perm_import) <- c("measure_name", "result")
  want_pull <- grep('.permutation.importance', perm_import$measure_name)
  final_perm_import <- perm_import[want_pull,]
  perm_remove <- which(final_perm_import$result == min(final_perm_import$result))
  if (length(perm_remove > 1)) {
    perm_remove <- perm_remove[1]
  }
  remove_this <- final_perm_import$measure_name[perm_remove]
  remove_this <- gsub('.permutation.importance', '', remove_this)

  # removing variable from analysis
  remove_env <- which(names(wanted_envs) %in% remove_this) 
  new_envs <- wanted_envs[[-remove_env]]
  
  return(new_envs)
}

# var_select function to iterate var_remove until all vifs are below 5
var_select <- function(wanted_envs, wanted_occs, wanted_bg) {
  
  final_env <- wanted_envs
  t <- max(usdm::vif(final_env)$VIF)
  
  while(t > 5) {
    final_env <- var_remove(wanted_envs = final_env, wanted_occs = wanted_occs,
                            wanted_bg = wanted_bg)
    if(dim(final_env)[3] < 3) {
      t <- 4
    } else {
      t <- max(usdm::vif(final_env)$VIF)
    }
  }
  
  return(final_env)
}
 
# SDM pipeline
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

tryCatch(
  {
  # reading in occurrences and accessible area
  occs <- read.csv(paste0('./Occurrences/', model_species[task_id], '.csv'))
  accessible_area <- st_read(dsn = './Acc_Area', layer = model_species[task_id])
  accessible_area <- st_set_crs(accessible_area, target_crs)
  
  # cropping environmental rasters
  envs_cropped <- crop(wanted_clim, accessible_area)
  envs_clipped <- mask(envs_cropped, accessible_area)
  
  # producing background
  acc_size <- ncell(envs_clipped[[1]])
  if (acc_size < 10000) {
    bgs <- bgs <- randomPoints(raster(envs_clipped[[1]]), n = round(acc_size*0.9))
  } else {
    bgs <- randomPoints(raster(envs_clipped[[1]]), n = 10000)
  }
  colnames(bgs) <- colnames(occs[,2:3])
  
  # variable selection
  envs_clipped_final <- var_select(wanted_envs = envs_clipped, wanted_occs = occs[,2:3],
                                   wanted_bg = bgs)
  
  # model training and validation
  wanted_kfold <- get.randomkfold(occs = occs[,2:3], bg = bgs, kfolds = 5)
  user.grp_kfold <- list(occs.grp = wanted_kfold$occs.grp, bg.grp = wanted_kfold$bg.grp)
  mod1 <- ENMevaluate(occs = occs[,2:3], envs = envs_clipped_final, bg = bgs,
                      algorithm = "maxnet",
                      tune.args = list(fc = c("L","LQ","LQH"), 
                                       rm = 1:3),
                      partitions = 'user', user.grp = user.grp_kfold)
  if (any(eval.results(mod1)$auc.val.avg) < 0.7) {
    select_this <- which(eval.results(mod1)$auc.val.avg == 
                           max(eval.results(mod1)$auc.val.avg, na.rm = T))
    if(length(select_this) > 1) {
      select_this <- select_this[1]
    }
    best_sdm <- mod1@predictions[[select_this]]
    best_sdm_tune <- mod1@tune.settings[select_this,]
  } else {
    select_this <- which(eval.results(mod1)$AICc == min(eval.results(mod1)$AICc, na.rm = T))
    if(length(select_this) > 1) {
      select_this <- select_this[1]
    }
    best_sdm <- mod1@predictions[[select_this]]
    best_sdm_tune <- mod1@tune.settings[select_this,]
  }
  write.table(c(best_sdm_tune, names(envs_clipped_final)), 
              file = paste0('./sdm_parameters/', model_species[task_id], '.txt'),
              row.names = F, col.names = F)
  
  # thresholding
  percentiles <- c(0.99, 0.95, 0.9)
  presences <- occs
  coordinates(presences) <- ~Longitude+Latitude
  proj4string(presences) <- target_crs
  absences <- data.frame(bgs)
  coordinates(absences) <- ~Longitude+Latitude
  proj4string(absences) <- target_crs
  sdm_suit <- extract(best_sdm, presences)[order(extract(best_sdm, presences), na.last = NA,
                                                 decreasing = T)]
  thresh_sdm_list <- vector("list", length = length(percentiles))
  tss_list <- vector(length = length(percentiles))
  for (j in seq_along(percentiles)) {
    thresh <- round(length(sdm_suit)*percentiles[j])
    thresh_val <- sdm_suit[thresh]
    thresholded_sdm <- reclassify(best_sdm, rcl = c(0, thresh_val, 0, thresh_val, 1, 1))
    thresh_sdm_list[[j]] <- thresholded_sdm
    real_abs <- extract(thresholded_sdm, absences)
    specificity <- length(which(real_abs == 0))/length(real_abs)
    sensitivity <- percentiles[j]
    tss <- (sensitivity + ((1/3)*specificity)) - 1
    tss_list[j] <- tss
  }
  
  choose_final <- which(tss_list == max(tss_list))
  final_sdm <- thresh_sdm_list[[choose_final]]
  names(final_sdm) <- model_species[task_id]
  
  # write raster to file
  writeRaster(final_sdm, file = paste0('./final_sdm_outputs/', model_species[task_id], '.tif'), 
              format = "GTiff", overwrite = T)

  }, 
  error=function(err){
    message('On iteration ',i, ' there was an error: ',err)
    file.copy(from = paste0("./Occurrences/", model_species[task_id], ".csv"), 
              to = paste0("./error_species/", model_species[task_id], ".csv"))
  }
)
