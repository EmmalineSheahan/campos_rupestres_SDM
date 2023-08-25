# Cleaning and combinging spatial occurrence datasets for the Poales group

library(sf)
library(terra)
library(dismo)
library(CoordinateCleaner)
library(dplyr)

# reading in datasets
poales_gbif <- read.csv('./data/poales_gbif.csv')
poales_idigbio <- read.csv('./data/poales_idigbio.csv')

# examining datasets
dim(poales_gbif)
names(poales_gbif)
poales_gbif$species

dim(poales_idigbio)
names(poales_idigbio)

# pulling wanted columns
poales_gbif <- poales_gbif %>% dplyr::select(species, decimalLatitude, decimalLongitude, 
                                      eventDate)
poales_sp_list_gbif <- unique(poales_gbif$species)

# fixing dates
poales_idigbio <- poales_idigbio %>% dplyr::select(two_word_species, 
                                                   data.dwc.decimalLatitude,
                                            data.dwc.decimalLongitude, data.dwc.month, 
                                            data.dwc.day, data.dwc.year)

# idigbio month
poales_idigbio$data.dwc.month <- gsub(".0", '', poales_idigbio$data.dwc.month)
for (i in 1:length(poales_idigbio$data.dwc.month)) {
  poales_idigbio$data.dwc.month[i] <- as.numeric(poales_idigbio$data.dwc.month[i])
  if (poales_idigbio$data.dwc.month[i] > 12 || is.na(poales_idigbio$data.dwc.month[i]) ||
      poales_idigbio$data.dwc.month[i] == 0) {
    poales_idigbio$data.dwc.month[i] <- NA
  } else if (poales_idigbio$data.dwc.month[i] == 1) {
    poales_idigbio$data.dwc.month[i] <- "01"
  } else if (poales_idigbio$data.dwc.month[i] == 2) {
    poales_idigbio$data.dwc.month[i] <- "02"
  } else if (poales_idigbio$data.dwc.month[i] == 3) {
    poales_idigbio$data.dwc.month[i] <- "03"
  } else if (poales_idigbio$data.dwc.month[i] == 4) {
    poales_idigbio$data.dwc.month[i] <- "04"
  } else if (poales_idigbio$data.dwc.month[i] == 5) {
    poales_idigbio$data.dwc.month[i] <- "05"
  } else if (poales_idigbio$data.dwc.month[i] == 6) {
    poales_idigbio$data.dwc.month[i] <- "06"
  }
}

#idigbio day
poales_idigbio$data.dwc.day <- gsub(".0", '', poales_idigbio$data.dwc.day)
for (i in 1:length(poales_idigbio$data.dwc.day)) {
  poales_idigbio$data.dwc.day[i] <- as.numeric(poales_idigbio$data.dwc.day[i])
  if (poales_idigbio$data.dwc.day[i] > 31 || is.na(poales_idigbio$data.dwc.day[i])) {
    poales_idigbio$data.dwc.day[i] <- NA
  }
}

# idigbio year
poales_idigbio$data.dwc.year <- gsub(".0", '', poales_idigbio$data.dwc.year)
poales_idigbio$date <- vector(length = length(poales_idigbio$data.dwc.month))
for (i in 1:length(poales_idigbio$data.dwc.month)) {
  if (is.na(poales_idigbio$data.dwc.month[i]) || is.na(poales_idigbio$data.dwc.day[i]) ||
      is.na(poales_idigbio$data.dwc.year[i])) {
        poales_idigbio$date[i] <- NA 
  } else {
        poales_idigbio$date[i] <- paste0(poales_idigbio$data.dwc.year[i], '-', 
                                   poales_idigbio$data.dwc.month[i], '-',
                                   poales_idigbio$data.dwc.day[i])
  }

}
