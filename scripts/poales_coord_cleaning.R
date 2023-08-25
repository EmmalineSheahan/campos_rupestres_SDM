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
library(spThin)

# creating land
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
land_poly <- as_Spatial(land)

# read in raw data
poales <- read.csv('./data/poales_data_campos_rupestres_2.csv')

# correcting synonyms and removing incorrect classifications
correct_syn <- read.csv('./data/poales_accepted_with_synomyns_minus_wrong_taxa.csv')
for (i in 1:nrow(correct_syn)) {
  t <- which(poales$two_word_species == correct_syn$synonyms[i])
  if (length(t) > 0) {
  poales$two_word_species[t] <- correct_syn$accepted_names[i]
  }
}

keep <- which(poales$two_word_species %in% correct_syn$accepted_names)
poales <- poales[keep,]

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

# creating list for nat ranges
nat_ranges <- list.files('./nat_ranges')
nat_ranges <- gsub(".dbf", '', nat_ranges)
nat_ranges <- gsub(".prj", '', nat_ranges)
nat_ranges <- gsub(".shp", '', nat_ranges)
nat_ranges <- gsub(".shx", '', nat_ranges)
nat_ranges <- unique(nat_ranges)

poales_names_list <- unique(poales_new$two_word_species)
poales_titles <- gsub("_", ' ', poales_names_list)

pdf('./plots/occ_clean_all_poales.pdf')
for (i in 1:length(poales_names_list)) {
  wanted_spec <- poales_new %>% filter(two_word_species == poales_names_list[i])
  wanted_spec2 <- wanted_spec %>% filter(!(is.na(geodeticDatum)))
  wanted_spec3 <- wanted_spec2 %>% filter(!(is.na(decimalLongitude)))
  wanted_spec4 <- wanted_spec3 %>% filter(!(is.na(decimalLatitude)))
  wanted_spec4 <- wanted_spec4[which(!(duplicated(wanted_spec4$Date))),]
  if (dim(wanted_spec4)[1] == 0) {
    print("No Coordinates")
  } else {
  sp_list <- vector("list", length = nrow(wanted_spec4))
  for (j in 1:nrow(wanted_spec4)) {
    test_sp <- wanted_spec4[j,]
    test_sp$decimalLongitude <- as.numeric(test_sp$decimalLongitude)
    test_sp$decimalLatitude <- as.numeric(test_sp$decimalLatitude)
    if (is.na(test_sp$geodeticDatum) & 
        test_sp$decimalLatitude > -180 & 
        test_sp$decimalLatitude < 180 & 
        test_sp$decimalLongitude > -90 &
        test_sp$decimalLongitude < 90) {
      test_sp$geodeticDatum <- "WGS84"
    } 
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
    } else if (test_sp$decimalLatitude > -180 & 
               test_sp$decimalLatitude < 180 & 
               test_sp$decimalLongitude > -90 &
               test_sp$decimalLongitude < 90) {
      coordinates(test_sp) <- ~decimalLongitude+decimalLatitude
      proj4string(test_sp) <- target_crs
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
  if (length(sp_all) < 1) {
    print("No Geodetic datum")
  } else if (length(sp_all) < 3) {
    print("only one record") 
  } else {
  sp_all <- data.frame(sp_all)
  sp_all_nodup <- unique(sp_all)
  species <- rep(poales_names_list[i], times = nrow(sp_all_nodup))
  sp_all_nodup <- cbind(species, sp_all_nodup)
  sp_all_nosea <- cc_sea(sp_all_nodup, lon = "decimalLongitude", lat = "decimalLatitude")
  if (poales_names_list[i] %in% nat_ranges) {
    want_range <- shapefile(paste0('./nat_ranges/', poales_names_list[i], '.shp'))
    sp_all_natrange <- cc_iucn(sp_all_nosea, want_range, lon = "decimalLongitude",
                             lat = "decimalLatitude", species = "species")
  } else {
    sp_all_natrange <- sp_all_nosea
  }
  if (nrow(sp_all_natrange) > 20) {
    sp_all_noout <- cc_outl(sp_all_natrange, lon = "decimalLongitude", 
                          lat = "decimalLatitude", species = "species", 
                          method = "distance", value = "clean", 
                          tdi = 1000)
  } else {
    sp_all_noout <- sp_all_natrange
  }
  if (nrow(sp_all_noout) == 0) {
    print("all outside natural range")
  } else {
  sp_thinned <- thin(loc.data = sp_all_noout, lat.col = "decimalLatitude", 
                     long.col = "decimalLongitude", spec.col = "species", thin.par = 10,
                     reps = 1, write.files = F, locs.thinned.list.return = T)
  sp_all_thinned <- sp_thinned[[1]]
  sp_occs <- sp_all_thinned
  coordinates(sp_occs) <- ~Longitude+Latitude
  plot(land_poly, col = NA, main = poales_titles[i],)
  plot(sp_occs, col = "red", pch = 19, cex = 0.25, add = T)
  if (poales_names_list[i] %in% nat_ranges) {
    plot(want_range, col = NA, border = "blue", add = T)
  }
  species <- rep(poales_names_list[i], times = nrow(sp_all_thinned))
  sp_all_thinned <- cbind(sp_all_thinned, species)
  write.csv(sp_all_thinned, file = paste0('./Occurrences/', poales_names_list[i], '.csv'))
  }
  }
  }
}
dev.off()

all_poales <- list.files('./Occurrences')
all_poales <- gsub('.csv', '', all_poales)
all_poales <- gsub('_', ' ', all_poales)
write.table(all_poales, file = './data/all_poales_species_list.txt', row.names = F,
            col.names = F)
