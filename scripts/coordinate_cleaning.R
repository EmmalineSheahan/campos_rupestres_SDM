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
poales_gbif <- poales_gbif %>% select(species, decimalLatitude, decimalLongitude, 
                                      eventDate)
poales_sp_list_gbif <- unique(poales_gbif$species)

# fixing dates
poales_idigbio <- poales_idigbio %>% select(two_word_species, data.dwc.decimalLatitude,
                                            data.dwc.decimalLongitude, data.dwc.month, 
                                            data.dwc.day, data.dwc.year)
poales_idigbio$data.dwc.month <- gsub(".0", '', poales_idigbio$data.dwc.month)
for (i in 1:length(poales_idigbio$data.dwc.month)) {
  if (as.numeric(poales_idigbio$data.dwc.month[i]) > 12) {
    poales_idigbio$data.dwc.month[i] <- NA
  }
}

poales_idigbio$data.dwc.day <- gsub(".0", '', poales_idigbio$data.dwc.day)
poales_idigbio$data.dwc.year <- gsub(".0", '', poales_idigbio$data.dwc.year)
poales_idigbio$date <- vector(length = length(poales_idigbio$data.dwc.month))
for (i in 1:length(poales_idigbio$data.dwc.month)) {
  poales_idigbio$date[i] <- paste0(poales_idigbio$data.dwc.year[i], '-', 
                                   poales_idigbio$data.dwc.month[i], '-',
                                   poales_idigbio$data.dwc.day[i])
}
