# coordinate cleaning poales

library(dplyr)
library(CoordinateCleaner)
library(terra)
library(sf)
library(dismo)
library(rgeos)
library(maptools)
library(rgdal)
library(rnaturalearth)

# creating land
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# read in raw data
poales <- read.csv('./data/poales_data_campos_rupestres.csv')

# investigate
names(poales)
length(unique(poales$two_word_species))
poales$year

# creating singular date column
poales$Date <- vector(length = length(poales$year))
for (i in 1:length(poales$year)) {
  poales$Date[i] <- paste0(poales$month[i], '-', poales$day[i], '-', poales$year[i])
}

# pulling needed columns
poales_new <- poales %>% dplyr::select(decimalLongitude, decimalLatitude, geodeticDatum,
                                       Date, two_word_species)

# adjusting two word name to avoid character string complications
poales_new$two_word_species <- gsub(' ', '_', poales_new$two_word_species)

# making sure all coordinates are numeric
for (i in 1:length(poales_new$decimalLatitude)) {
  poales_new$decimalLatitude[i] <- as.numeric(poales_new$decimalLatitude[i])
  poales_new$decimalLongitude[i] <- as.numeric(poales_new$decimalLongitude[i])
}

# creating needed crs' for transforming
epsg <- make_EPSG()
target_crs <- epsg %>% filter(code == 4326)
target_crs <- target_crs$prj4

nad27_crs <- epsg %>% filter(code == 4267)
nad27_crs <- nad27_crs$prj4

nad83_crs <- epsg %>% filter(code == 4269)
nad83_crs <- nad83_crs$prj4

merc_crs <- epsg %>% filter(code == 3857)
merc_crs <- merc_crs$prj4

# creating and cleaning spatial points dataframes per species
dir.create('./Occurrences')

poales_names_list <- unique(poales_new$two_word_species)
for (i in 1:length(poales_names_list)) {
  wanted_spec <- poales_new %>% filter(two_word_species == poales_names_list[i])
  wanted_spec2 <- wanted_spec %>% filter(!(is.na(geodeticDatum)))
  wanted_spec3 <- wanted_spec2 %>% filter(!(is.na(decimalLongitude)))
  wanted_spec4 <- wanted_spec3 %>% filter(!(is.na(decimalLatitude)))
  sp_list <- vector("list", length = nrow(wanted_spec4))
  for (j in 1:nrow(wanted_spec4)) {
    test_sp <- wanted_spec4[j,]
    test_sp$decimalLongitude <- as.numeric(test_sp$decimalLongitude)
    test_sp$decimalLatitude <- as.numeric(test_sp$decimalLatitude)
    if (test_sp$geodeticDatum == "WGS84" || 
        test_sp$geodeticDatum == "World Geodetic System 1984" ||
        test_sp$geodeticDatum == "EPSG:4326" ||
        test_sp$geodeticDatum == "WGS-84" ||
        test_sp$geodeticDatum == "WGS 84" ||
        test_sp$geodeticDatum == "WGS 1984" ||
        test_sp$geodeticDatum == "WGs84" ||
        test_sp$geodeticDatum == "wgs84") {
      coordinates(test_sp) <- ~decimalLongitude+decimalLatitude
      proj4string(test_sp) <- target_crs
    } else if (test_sp$geodeticDatum == "NAD27" ||
               test_sp$geodeticDatum == "NAD 1927" ||
               test_sp$geodeticDatum == "NAD 27") {
      coordinates(test_sp) <- ~decimalLongitude+decimalLatitude
      proj4string(test_sp) <- nad27_crs
      test_sp <- spTransform(test_sp, CRSobj = target_crs)
    } else if (test_sp$geodeticDatum == "NAD 83" ||
               test_sp$geodeticDatum == "NAD83" ||
               test_sp$geodeticDatum == "Nad83" ||
               test_sp$geodeticDatum == "Nad 83") {
      coordinates(test_sp) <- ~decimalLongitude+decimalLatitude
      proj4string(test_sp) <- nad83_crs
      test_sp <- spTransform(test_sp, CRSobj = target_crs)
    } else if (test_sp$geodeticDatum == "Google Maps USGS Maps Google Earth") {
      coordinates(test_sp) <- ~decimalLongitude+decimalLatitude
      proj4string(test_sp) <- merc_crs
      test_sp <- spTransform(test_sp, CRSobj = target_crs)
    } else {
      test_sp <- NA
    }
    if (is.na(test_sp)) {
      sp_list[[j]] <- NA
    } else {
    sp_list[[j]] <- test_sp@coords
    }
  }
  sp_all <- do.call(rbind, sp_list)
  sp_all <- sp_all[complete.cases(sp_all),]
}
