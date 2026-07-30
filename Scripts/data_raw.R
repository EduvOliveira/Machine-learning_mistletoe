#' ---
# Disentangling biotic and abiotic drivers of a neotropical mistletoe 
# By Eduardo V. S. Oliveira
# 28/07/2026
#' ---

## Data acquisition and organization

library(rgbif)
library(plyr)
library(dplyr)
library(raster)
library(BIEN)
library(dismo)
library(sp)
library(sf)
library(flora)
library(terra)
library(stringr)
library(plantR)

# Data directly from GBIF 

avian.spp<-occ_search(scientificName=c("Tersina viridis viridis", "Elaenia cristata","Thraupis sayaca sayaca","Schistochlamys ruficapillus ruficapillus"),hasCoordinate = TRUE, fields = c("scientificName","decimalLatitude","decimalLongitude","country"))

avian.spp$"Tersina viridis"$data

bdspp<-rbind.fill(avian.spp$"Tersina viridis viridis"$data,avian.spp$"Elaenia cristata"$data,
                  avian.spp$"Thraupis sayaca sayaca"$data,avian.spp$"Schistochlamys ruficapillus ruficapillus"$data)

setwd("J:/path_to_the_folder")

occ_aves <- rename(bdspp, spp = scientificName, lon = decimalLongitude, lat = decimalLatitude)

occ_birds<-occ_aves[,1:3]

# Plant occurrence data from BIEN

species_vector<-c("Psittacanthus robustus","Salvertia convallariodora", "Qualea dichotoma", "Qualea grandiflora", "Qualea multiflora", "Qualea parviflora", "Vochysia cinnamomea", "Vochysia elliptica", "Vochysia rufa", "Vochysia tucanorum", "Vochysia thyrsoidea","Trembleya laniflora","Qualea cordata","Miconia ferruginata","Miconia albicans","Diplusodon hirsutus","Coccoloba cereifera","Campomanesia adamantium")

get.taxa(species_vector)

occ_data<-BIEN_occurrence_species(species_vector) 

occ_plants<-occ_data[,1:3]

occ_plants<- rename(occ_plants, spp = scrubbed_species_binomial, lat = latitude, lon = longitude)

occ_veg<-occ_plants %>% relocate(lon, .before = lat)

# Plant occurrence data from GBIF

plant.spp<-occ_search(scientificName=c("Psittacanthus robustus","Salvertia convallariodora", "Qualea dichotoma", "Qualea grandiflora", "Qualea multiflora", "Qualea parviflora", "Vochysia cinnamomea", "Vochysia elliptica", "Vochysia rufa", "Vochysia tucanorum", "Vochysia thyrsoidea","Trembleya laniflora","Qualea cordata","Miconia ferruginata","Miconia albicans","Diplusodon hirsutus","Coccoloba cereifera","Campomanesia adamantium"),hasCoordinate = TRUE, fields = c("scientificName","decimalLatitude","decimalLongitude","country"))

occ.spp<-rbind.fill(plant.spp$"Psittacanthus robustus"$data,plant.spp$"Salvertia convallariodora"$data,plant.spp$"Qualea dichotoma"$data,plant.spp$"Qualea grandiflora"$data,plant.spp$"Qualea multiflora"$data,plant.spp$"Qualea parviflora"$data,plant.spp$"Vochysia cinnamomea"$data,plant.spp$"Vochysia elliptica"$data,plant.spp$"Vochysia rufa"$data,plant.spp$"Vochysia tucanorum"$data,plant.spp$"Vochysia thyrsoidea"$data,plant.spp$"Trembleya laniflora"$data,plant.spp$"Qualea cordata"$data,plant.spp$"Miconia ferruginata"$data,plant.spp$"Miconia albicans"$data,plant.spp$"Diplusodon hirsutus"$data,plant.spp$"Coccoloba cereifera"$data,plant.spp$"Campomanesia adamantium"$data)

occ_angio <- rename(occ.spp, spp = scientificName, lon = decimalLongitude, lat = decimalLatitude)

occ_angio<-occ_angio[,1:3]

# Additionally, obtaining plant and bird records from SpeciesLink

occ_splink <- lapply(lista, function(lista) rspeciesLink(species = lista,Coordinates = "Yes"))

# Name the objects in the list using the species name

names(occ_splink) = lista

# Convert the list into a dataframe

occs_splink <- bind_rows(occs_splink, .id = "original_search")

# Putting it all together

bdtotal<-rbind.fill(occ_birds, occ_veg,occ_angio,occ_splink)

# Standardizing the names

occ <- bdtotal %>%
  mutate(
    species_original = spp,
    
        spp = str_replace_all(spp, "cf\\.|aff\\.|sp\\.|spp\\.", ""),spp = word(spp, 1, 2),
    spp = str_squish(spp)
  )

write.csv(occ, "occ_spp_raw.csv")

# plotting the points with the South America polygon

Am_sul<-shapefile(choose.files())
plot(Am_sul)
points(occ$lon, occ$lat, col='red', pch=20, cex=0.75) 

# Removing duplicate records

occ<-occ[,1:3]

dups <- duplicated(occ[, c('spp', 'lon', 'lat')]) sum(dups) 
sp_sub <- occ[!dups, ] 
dim(sp_sub)

plot(Am_sul)
points(sp_sub$lon, sp_sub$lat, col='red', pch=20, cex=0.75) 

write.csv(sp_sub, "occ_sp_sub.csv")

# Removing records outside of South America

occ <- read.csv("occ_sp_sub.csv", stringsAsFactors = FALSE)

occ<-occ[,2:4]

# Standardizing the names
occ <- occ %>%
  rename(
    species = spp,
    lon = lon,
    lat = lat
  ) %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(!is.na(lon), !is.na(lat))


# Convert to spatial object
occ_sf <- st_as_sf(
  occ,
  coords = c("lon", "lat"),
  crs = 4326 # WGS84
)

plot(occ_sf)

# Load South America shapefile

dir()
sa <- st_read("South_America.shp")

sa <- st_union(sa)

plot(sa)

sa <- st_transform(sa, crs = 4326)

# Filter points inside the polygon
occ_in_sa <- occ_sf[st_within(occ_sf, sa, sparse = FALSE), ]

plot(occ_in_sa)

# Remove geometry and return to dataframe

occ_final <- occ_in_sa %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# Save
write.csv(occ_final, "occ_south_america.csv", row.names = FALSE)

# Mapping the data again

ocu <- read.csv("occ_south_america.csv", stringsAsFactors = FALSE)

plot(Am_sul)
points(ocu$lon, ocu$lat, col='red', pch=20, cex=0.75)


# Addressing sampling biases
# Select only one record per 1-degree cell

# Create a 1-degree grid
bb <- st_bbox(occ_in_sa)

grid <- st_make_grid(
  st_as_sfc(bb),
  cellsize = 1,     
  square = TRUE
)

grid <- st_sf(grid_id = 1:length(grid), geometry = grid)

plot(grid)

# Point-grid intersection
occ_grid <- st_join(occ_in_sa, grid, join = st_within)

# Thinning
occ_thinned <- occ_grid %>%
  group_by(species, grid_id) %>%
  slice(1) %>%   
  ungroup()

occ_thinned <- occ_thinned %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  )

plot(occ_thinned)

# Remove geometry

occ_thin <- occ_thinned %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# Save
write.csv(occ_thin, "occ_thinned_1deg.csv", row.names = FALSE)

# Organize the data

levels(as.factor(occ_thin$species))

species_choice<-c("Psittacanthus robustus","Salvertia convallariodora", "Qualea dichotoma", "Qualea grandiflora", "Qualea multiflora", "Qualea parviflora", "Vochysia cinnamomea", "Vochysia elliptica", "Vochysia rufa", "Vochysia tucanorum", "Vochysia thyrsoidea","Trembleya laniflora","Qualea cordata","Miconia ferruginata","Miconia albicans","Diplusodon hirsutus","Coccoloba cereifera","Campomanesia adamantium","Tersina viridis", "Elaenia cristata","Thraupis sayaca","Schistochlamys ruficapillus")


occ_grid <- occ_thin %>%
  filter(species %in% species_choice)

levels(as.factor(occ_grid$species))

occ_grid$n<-rep(1,nrow(occ_grid))

names(occ_grid)

occ_grid<-occ_grid %>% dplyr::select (species, lon, lat, n)

write.table(occ_grid, "occ_clean.txt")

##THE END##

rm(list=ls())
