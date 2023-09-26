# creating PAMs

library(dplyr)
library(sf)
library(raster)
library(terra)
library(ape)
library(rnaturalearth)
library(reshape2)

# all species list
# I would copy the rasters from the species we made hulls for into the final_sdm_outputs
# directory so everything is in one place

all_species <- list.files('./final_sdm_outputs')
all_species <- gsub('.tif', '', all_species)

# read in the tree
full_tree <- read.nexus('./data/pruned_rooted_poales_for_PD.tre')

# prune the tree to the list of species we have rasters for
tree_names <- full_tree$tip.label

# you'll need to make sure the tree_names are in the same format (Genus_species) as species
# are listed in all_species, you can use the function gsub to alter that
clean_names <- gsub(' ', '_', tree_names)
full_tree$tip.label <- clean_names
remove_these <- which(!(full_tree$tip.label %in% all_species))
if(length(remove_these) != 0) {
  final_tree <- drop.tip(full_tree, remove_these)
} else {
  final_tree <- full_tree
}

# read in species stack, change to raster stack
species_stack_1 <- rast('./data/all_species_stack.tif')

species_stack_2.5 <- aggregate(species_stack_1, fact = 2.5)

species_stack_5 <- aggregate(species_stack_1, fact = 5)

species_stack_list <- vector("list", length = dim(species_stack_1)[3])
for (i in 1:dim(species_stack_1)[3]) {
  species_stack_list[[i]] <- raster(species_stack_1[[i]])
}
cr_stack_1 <- stack(species_stack_list)

species_stack_list <- vector("list", length = dim(species_stack_2.5)[3])
for (i in 1:dim(species_stack_2.5)[3]) {
  species_stack_list[[i]] <- raster(species_stack_2.5[[i]])
}
cr_stack_2.5 <- stack(species_stack_list)

species_stack_list <- vector("list", length = dim(species_stack_5)[3])
for (i in 1:dim(species_stack_5)[3]) {
  species_stack_list[[i]] <- raster(species_stack_5[[i]])
}
cr_stack_5 <- stack(species_stack_list)

# creating row remove of false NA's
# read in richness raster for base
cr_richness_1 <- raster('./data/campos_rupestres_richness.tif')
row_remove_1 <- which(is.na(values(cr_richness_1)))

cr_richness_2.5 <- aggregate(cr_richness_1, fact = 2.5)
row_remove_2.5 <- which(is.na(values(cr_richness_2.5)))

cr_richness_5 <- aggregate(cr_richness_1, fact = 5)
row_remove_5 <- which(is.na(values(cr_richness_5)))

# presence absence matrix function
pam_maker <- function(wanted_stack) {
  mat <- matrix(NA, ncol = dim(wanted_stack)[3], 
                nrow = ncell(wanted_stack[[1]]) + 1)
  for (i in 1:dim(wanted_stack)[3]) {
    sname <- names(wanted_stack[[i]])
    mat[1,i] <- sname
    mat[2:nrow(mat),i] <- as.numeric(values(wanted_stack[[i]]))
  }
  mat <- data.frame(mat)
  names(mat) <- lapply(mat[1, ], as.character)
  mat <- mat[-1,]
  return(mat)
}

cr_pam_1 <- pam_maker(cr_stack_1)
cr_pam_2.5 <- pam_maker(cr_stack_2.5)
cr_pam_5 <- pam_maker(cr_stack_5)

# coordinates matrix function
coord_maker <- function(wanted_stack) {
  mat <- matrix(NA, ncol = 2, 
                nrow = ncell(wanted_stack[[1]]) + 1)
  lay <- wanted_stack[[1]]
  rclmat <- matrix(c(NA, 0), ncol = 2, byrow = T)
  lay <- reclassify(lay, rclmat)
  pts <- rasterToPoints(lay, spatial = T)
  mat[2:nrow(mat),1] <- as.numeric(pts@coords[,1])
  mat[2:nrow(mat),2] <- as.numeric(pts@coords[,2])
  mat[1,1] <- "lon"
  mat[1,2] <- "lat"
  mat <- data.frame(mat)
  names(mat) <- lapply(mat[1, ], as.character)
  mat <- mat[-1,]
  return(mat)
}

cr_coords_1 <- coord_maker(cr_stack_1)
cr_coords_2.5 <- coord_maker(cr_stack_2.5)
cr_coords_5 <- coord_maker(cr_stack_5)

# matching pam to tree with pam_matcher function
# if you want absences included in the pam, set include_absence = T
# biodiverse wants pams with just presences
pam_matcher <- function(wanted_pam, wanted_tree, 
                        wanted_coords, wanted_row_remove, include_absence = F) {
  wanted_pam1 <- wanted_pam %>% dplyr::select(wanted_tree$tip.label)
  
  for (i in 1:ncol(wanted_pam1)) {
    wanted_pam1[,i] <- as.numeric(wanted_pam1[,i])
  }
  
  for (i in 1:ncol(wanted_coords)) {
    wanted_coords[,i] <- as.numeric(wanted_coords[,i])
  }
  
  wanted_pam1clean <- wanted_pam1[-wanted_row_remove,]
  wanted_coordsclean <- wanted_coords[-wanted_row_remove,]
  
  wanted_pam2 <- wanted_pam1clean %>% replace(is.na(.), 0)
  wanted_pam_new <- cbind(wanted_coordsclean, wanted_pam2)
  
  spsumsrow <- rowSums(wanted_pam_new[,3:ncol(wanted_pam_new)])
  droprow <- which(spsumsrow == 0)
  if (length(droprow) > 0) {
    wanted_pam_final1 <- wanted_pam_new[-droprow,]
  } else {
    wanted_pam_final1 <- wanted_pam_new
  }
  
  spsumscol <- colSums(wanted_pam_final1)
  dropcol <- which(spsumscol == 0)
  if (length(dropcol) > 0) {
    wanted_pam_final <- wanted_pam_final1[,-dropcol]
  } else {
    wanted_pam_final <- wanted_pam_final1
  }
  wanted_pam_final_melt <- melt(wanted_pam_final, id = c("lon", "lat"))
  if(include_absence == F) {
    wanted_pam_final_melt <- wanted_pam_final_melt %>% filter(value == 1)
  }
  return(wanted_pam_final_melt)
}

cr_pam_final_1 <- pam_matcher(cr_pam_1, final_tree, cr_coords_1, row_remove_1, 
                              include_absence = F)
write.csv(cr_pam_final_1, file = './data/cr_pam_final_1.csv')

cr_pam_final_2.5 <- pam_matcher(cr_pam_2.5, final_tree, cr_coords_2.5, row_remove_2.5, 
                              include_absence = F)
write.csv(cr_pam_final_2.5, file = './data/cr_pam_final_2.5.csv')

cr_pam_final_5 <- pam_matcher(cr_pam_5, final_tree, cr_coords_5, row_remove_5, 
                              include_absence = F)
write.csv(cr_pam_final_5, file = './data/cr_pam_final_5.csv')