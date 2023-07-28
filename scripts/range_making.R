# script to create range maps for the natural range cleaning function of CoordinateCleaner

library(sf)
library(sp)
library(rnaturalearth)
library(dplyr)
library(geobr)
library(rgeos)
library(terra)
library(rgdal)
library(raster)

# reading in data
incorrect_syn <- read.csv('./data/incorrect_syn_wcvp.csv')
wcvp_dist <- read.csv('./data/wcvp_taxon_poales_cr_dist_filtered_2.csv')
brazil_dist <- read.csv('./data/poales_distribution_brazil.csv')

# correcting names to match with poales data
wrong_names <- incorrect_syn$Accepted.by.WCVP.but.should.be.considered.synonyms
correct_names <- incorrect_syn$Synonyms.by.WCVP.but.should.be.considered.accepted

for (i in 1:length(wcvp_dist$scientfiicname)) {
  if (wcvp_dist$scientfiicname[i] %in% wrong_names) {
    wanted_num <- which(wcvp_dist$scientfiicname[i] %in% wrong_names)
    wcvp_dist$scientfiicname[i] <- correct_names[wanted_num]
  }
  wcvp_dist$scientfiicname[i] <- gsub(' ', '_', wcvp_dist$scientfiicname[i])
}

# creating directory for range maps
dir.create('./nat_ranges')
dir.create('./brazil_ranges')

# creating phytogeographic regions
biomes <- read_biomes(year = 2019, simplified = TRUE)
biomes_no_sea <- biomes[-7,]

biomes_st <- st_as_sf((biomes_no_sea), remove = FALSE, crs = 4326, agr = "constant")[1]
biomes_st$name_biome <- c("Amazonia", "Caatinga", "Cerrado", "Mata_Atlantica", "Pampa",
                          "Pantanal")

# creating crs for transformation
epsg <- make_EPSG()
target_crs <- epsg %>% filter(code == 4326)
target_crs <- target_crs$prj4

# adjusting locality to be readable for rnaturalearth
locality <- unique(wcvp_dist$locality)
world <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
country_names <- unique(world$name)
brazil <- ne_states(country = "Brazil", returnclass = "sf")

wcvp_dist$locality <- gsub(' Northeast', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' Southeast', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' Northwest', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' Southwest', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' North', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' South', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' East', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' West', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' West-Central', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' South-Central', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' East-Central', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' North-Central', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub(' Central', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('-Central', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Western ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Northern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Eastern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Southern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Central ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('North ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('East ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('South ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('West ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Northeastern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Northwestern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Southwestern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Southeastern ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Northeast ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Northwest ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Southwest ', '', wcvp_dist$locality)
wcvp_dist$locality <- gsub('Southeast ', '', wcvp_dist$locality)
locality <- unique(wcvp_dist$locality)

flag <- integer()
for (i in 1:length(locality)) {
  tryCatch(
    {
      testing <- ne_countries(country = locality[i], returnclass = "sf", scale = 50)
    }, 
    error=function(err){
      message('On iteration ',i, ' there was an error: ',err)
      flag <<-c(flag,i)
    }
  )
}

problem_loc <- locality[flag]
working_countries <- locality[-flag]
problem_loc[120] <- "Rhode Island"

statesflag <- integer()
for (i in 1:length(problem_loc)) {
  tryCatch(
    {
      testing <- ne_states(country = "United States of America", returnclass = "sf")[9]
      newtest <- testing %>% filter(name == problem_loc[i])
      plot(newtest)
    }, 
    error=function(err){
      message('On iteration ',i, ' there was an error: ',err)
      statesflag <<-c(statesflag,i)
    }
  )
}

working_us_states <- problem_loc[-statesflag]
still_problems <- problem_loc[statesflag]

china <- ne_states(country = "China", returnclass = "sf")[9]

chinaflag <- integer()
for (i in 1:length(still_problems)) {
  tryCatch(
    {
      testing <- ne_states(country = "China", returnclass = "sf")[9]
      newtest <- testing %>% filter(name == still_problems[i])
      plot(newtest)
    }, 
    error=function(err){
      message('On iteration ',i, ' there was an error: ',err)
      chinaflag <<-c(chinaflag,i)
    }
  )
}

working_chinese_states <- still_problems[-chinaflag]
still_problems <- still_problems[chinaflag]
still_problems[16] <- "New South Wales"

ausflag <- integer()
for (i in 1:length(still_problems)) {
  tryCatch(
    {
      testing <- ne_states(country = "Australia", returnclass = "sf")[9]
      newtest <- testing %>% filter(name == still_problems[i])
      plot(newtest)
    }, 
    error=function(err){
      message('On iteration ',i, ' there was an error: ',err)
      ausflag <<-c(ausflag,i)
    }
  )
}

working_aus_states <- still_problems[-ausflag]
still_problems <- still_problems[ausflag]

canflag <- integer()
for (i in 1:length(still_problems)) {
  tryCatch(
    {
      testing <- ne_states(country = "Canada", returnclass = "sf")[9]
      newtest <- testing %>% filter(name == still_problems[i])
      plot(newtest)
    }, 
    error=function(err){
      message('On iteration ',i, ' there was an error: ',err)
      canflag <<-c(canflag,i)
    }
  )
}

working_can_states <- still_problems[-canflag]
still_problems <- still_problems[canflag]

safflag <- integer()
for (i in 1:length(still_problems)) {
  tryCatch(
    {
      testing <- ne_states(country = "South Africa", returnclass = "sf")[9]
      newtest <- testing %>% filter(name == still_problems[i])
      plot(newtest)
    }, 
    error=function(err){
      message('On iteration ',i, ' there was an error: ',err)
      safflag <<-c(safflag,i)
    }
  )
}

working_saf_states <- still_problems[-safflag]
still_problems <- still_problems[safflag]

change_to_problems <- c(NA, "Panama", "Antigua and Barbuda", "Mexico", 
                        "Trinidad and Tobago", NA, "Burkina Faso", NA, 
                        "Central African Republic", "Democratic Republic of the Congo",
                        "Guinea Bissau", NA, "Indonesia", NA, "Papua New Guinea", NA, 
                        NA, NA, NA, NA, "United Republic of Tanzania", NA, 
                        "Democratic Republic of the Congo", "India", NA, 
                        "Wallis and Futuna",
                        NA, "Namibia", "The Bahamas", "Cayman Islands", NA, NA, NA, 
                        "Ecuador", NA, NA, "Portugal", "Lithuania", "Russia", NA, 
                        "Czech Republic", "Greece", "Russia", "Kyrgyzstan", "Greece",
                        "Ukraine", "Lebanon", "Portugal", "Georgia", "Russia", "Italy",
                        "Italy", "Russia", "Tajikistan", "Georgia", "Turkey", "India",
                      "Montenegro", "Australia", "Indonesia", "Kiribati", "Australia",
                      "Costa Rica", "New Zealand", NA, "Japan", "South Korea", 
                      "Indonesia", "Malaysia", "Indonesia", NA, "Marshall Islands", 
                      "Japan", "India", NA, "Japan", "Indonesia", "Indonesia", NA,
                      "India", NA, "Papua New Guinea", "Saudi Arabia", NA, NA, 
                      "Mozambique", "kiribati", NA, "Yemen", "Solomon Islands", "China",
                      "Turks and Caicos Islands", NA, "Marshall Islands", NA, NA, "China",
                      "Rhode Island", "New South Wales")

still_problems <- c(still_problems, "Rhode I.", "New Wales")

for (i in 1:length(wcvp_dist$locality)) {
  if(wcvp_dist$locality[i] %in% still_problems) {
    change_this <- which(still_problems %in% wcvp_dist$locality[i])
    wcvp_dist$locality[i] <- change_to_problems[change_this]
  }
}

# Separating list into Brazil species and non Brazil species
brazil_species <- brazil_dist$TAXON
brazil_species <- gsub(' ', '_', brazil_species)
for (i in 1:length(wcvp_dist$scientfiicname)) {
  loc <- wcvp_dist[i,] %>% dplyr::select(locality)
  if(length(loc) == 1 & loc == "Brazil" | is.na(loc)) {
    wcvp_dist[i,]$scientfiicname <- NA
  }
}
drop_spec <- which(is.na(wcvp_dist$scientfiicname))
wcvp_dist <- wcvp_dist[-drop_spec,]
wcvp_species <- wcvp_dist$scientfiicname
wcvp_species <- unique(wcvp_species)

# creating species range polygons for coordinate cleaning
for (i in seq_along(brazil_species)) {
  brazil_regions <- brazil_dist[, 4:9]
  reg_num <- sum(brazil_regions[i,])
  if (!(is.na(reg_num))) {
  region_names <- names(brazil_regions)
  sp_poly <- vector("list", length = length(region_names))
  for (j in seq_along(region_names)) {
    wanted_spec <- brazil_regions[i,]
    if (wanted_spec %>% dplyr::select(region_names[j]) == 1) {
      reg_poly <- biomes_st %>% filter(name_biome == region_names[j])
      reg_poly <- as_Spatial(reg_poly)
      reg_poly <- spTransform(reg_poly, CRSobj = target_crs)
      sp_poly[[j]] <- reg_poly
    } else {
      sp_poly[[j]] <- NA
    }
  }
  drop_na <- which(is.na(sp_poly))
  if (length(drop_na) > 0) {
    sp_poly_new <- sp_poly[-drop_na]
  } else {
    sp_poly_new <- sp_poly
  }
  if (length(sp_poly_new) == 1) {
    sp_poly_brazil <- sp_poly_new[[1]]
  } else {
    sp_poly_brazil <- do.call(raster::bind, sp_poly_new)
  }
  writeOGR(sp_poly_brazil, dsn = './brazil_ranges/', layer = brazil_species[i],
           driver = 'ESRI Shapefile', overwrite = T)
  if (brazil_species[i] %in% wcvp_species) {
    wcvp_wantspec <- wcvp_dist %>% filter(scientfiicname == brazil_species[i])
    wantspec_locality <- wcvp_wantspec %>% dplyr::select(locality)
    wantspec_locality <- unique(wantspec_locality$locality)
    add_poly <- vector("list", length = length(wantspec_locality))
     for (j in seq_along(wantspec_locality)) {
      if (wantspec_locality[j] %in% working_aus_states) {
        testing <- ne_states(country = "Australia", returnclass = "sf")[9]
        newtest <- testing %>% filter(name == wantspec_locality[j])
        newtest <- as_Spatial(newtest)
        newtest <- spTransform(newtest, CRSobj = target_crs)
        add_poly[[j]] <- newtest
      } else if (wantspec_locality[j] %in% working_chinese_states) {
        testing <- ne_states(country = "China", returnclass = "sf")[9]
        newtest <- testing %>% filter(name == wantspec_locality[j])
        newtest <- as_Spatial(newtest)
        newtest <- spTransform(newtest, CRSobj = target_crs)
        add_poly[[j]] <- newtest
      } else if (wantspec_locality[j] %in% working_can_states) {
        testing <- ne_states(country = "Canada", returnclass = "sf")[9]
        newtest <- testing %>% filter(name == wantspec_locality[j])
        newtest <- as_Spatial(newtest)
        newtest <- spTransform(newtest, CRSobj = target_crs)
        add_poly[[j]] <- newtest
      } else if (wantspec_locality[j] %in% working_saf_states) {
        testing <- ne_states(country = "South Africa", returnclass = "sf")[9]
        newtest <- testing %>% filter(name == wantspec_locality[j])
        newtest <- as_Spatial(newtest)
        newtest <- spTransform(newtest, CRSobj = target_crs)
        add_poly[[j]] <- newtest
      } else if (wantspec_locality[j] %in% working_us_states) {
        testing <- ne_states(country = "United states of America", returnclass = "sf")[9]
        newtest <- testing %>% filter(name == wantspec_locality[j])
        newtest <- as_Spatial(newtest)
        newtest <- spTransform(newtest, CRSobj = target_crs)
        add_poly[[j]] <- newtest
      } else {
        testing <- ne_countries(country = wantspec_locality[j], returnclass = "sf", 
                                scale = 50)[18]
        newtest <- as_Spatial(testing)
        newtest <- spTransform(newtest, CRSobj = target_crs)
        add_poly[[j]] <- newtest
      }
    }
  if (length(add_poly) == 1) {
    add_poly_full <- add_poly[[1]]
  } else {
    add_poly_full <- do.call(raster::bind, add_poly)
  }
  final_poly <- raster::bind(add_poly_full, sp_poly_brazil)
  final_poly@data$species <- rep(brazil_species[i], times = dim(final_poly@data)[1])
  writeOGR(final_poly, dsn = './nat_ranges/', layer = brazil_species[i],
           driver = 'ESRI Shapefile', overwrite = T)
  } else {
    sp_poly_brazil@data$species <- rep(brazil_species[i], 
                                       times = dim(sp_poly_brazil@data)[1])
    writeOGR(sp_poly_brazil, dsn = './nat_ranges/', layer = brazil_species[i],
             driver = 'ESRI Shapefile', overwrite = T)
  }
  }
}

# creating plots of natural ranges
nat_ranges <- list.files('./nat_ranges')
nat_ranges <- gsub(".dbf", '', nat_ranges)
nat_ranges <- gsub(".prj", '', nat_ranges)
nat_ranges <- gsub(".shp", '', nat_ranges)
nat_ranges <- gsub(".shx", '', nat_ranges)
nat_ranges <- unique(nat_ranges)

nat_titles <- gsub("_", ' ', nat_ranges)

pdf('./plots/natural_distribution_polygons.pdf')
for (i in 1:length(nat_ranges)) {
  range_poly <- shapefile(paste0('./nat_ranges/', nat_ranges[i], '.shp'))
  plot(land_poly, col = NA, main = nat_titles[i])
  plot(range_poly, col = "yellow", add = T)
}
dev.off()