#' ---
# Disentangling biotic and abiotic drivers of a neotropical mistletoe 
# By Eduardo V. S. Oliveira
# 28/01/2026
#' ---

install.packages(c('sf', 'terra', 'tidyverse','here', 'raster','openxlsx'))

# Load libraries

library(sf)
library(terra)
library(dplyr)
library(here)
library(raster)
library(openxlsx)

# Load data

here::i_am("grids_ML.R") 

occ<-read.csv(here("Data","occ_spp_raw.csv"))

head(occ)

p<-occ[,2:4]

head(p)

pontos_sf <- st_as_sf(p, coords = c("lon", "lat"), crs = 4326)

# Create buffers around the points and collapse them
buffer_dist <- 25000  # Buffer distance (25 km)

pontos_buffer <- pontos_sf %>%
  group_by(spp) %>%
  reframe(geometry = st_union(st_buffer(geometry, dist = buffer_dist))) %>%
  st_as_sf() 

pontos_buffer <- pontos_buffer %>%
  mutate(geometry = st_make_valid(geometry))

plot(pontos_buffer)

# Save the buffer with shapefile
st_write(pontos_buffer, here("data", "buffers_colapsados.shp"))

# Create a grid for South America

# Define the South America extent
ext <- ext(-85, -30, -60, 15)  

# Create the raster with the desired resolution (in degrees) and assign a CRS
resolucao <- 1  
grid <- rast(ext, resolution = resolucao, crs = "EPSG:4326")

# Convert the raster to an sf object
grid_sf <- as.polygons(grid) %>% st_as_sf()

# Intersection between collapsed buffers and the grid

# Filter valid buffers
buffers_validos <- buffers_colapsados %>% filter(!st_is_empty(geometry))

# Perform the intersection between the grid and the buffers
incidencia <- st_intersects(grid_sf, buffers_validos, sparse = FALSE)

# Create a presence/absence matrix
presenca_ausencia <- ifelse(incidencia, 1, 0)

# Naming the matrix columns based on the species
colnames(presenca_ausencia) <- paste0("_", unique(buffers_validos$spp))

# Inspect the presence/absence matrix
head(presenca_ausencia)

# Add the incidence matrix to the grid
grid_sf <- grid_sf %>%
  bind_cols(as.data.frame(presenca_ausencia))

# Save the grid with the presence/absence matrix as a shapefile
st_write(grid_sf, "grid_inc4.shp")

write.csv(presenca_ausencia, "mpa.csv")

# Drawing lots for the pseudo-absences

g<-shapefile("grid_inc4.shp")

a<-shapefile("South_America.shp")

plot(g)
plot(a, add=TRUE)

gi<-intersect(g,a)

plot(gi)

shapefile(gi,"grid_south_a.shp")

g<-st_read("grid_south_a.shp")

mpa<-dplyr::select(g, X_EcP.18, X_Psttcr,X_Qldcht,X_Qlgrnd,X_Qlmltf,X_Qlprvf,X_Slvrtc,X_Schsrr ,X_Trsnvv,X_Thrpss,X_Vchysc,X_Vchyse,X_Vchysr,X_Vchysth,X_Vchystc)

head(mpa)

# Randomly selecting the presences

p<-mpa%>%dplyr::filter(X_Psttcr == 1)

dim(p)

st_write(p, "grid_p.shp")

a<-mpa%>%dplyr::filter(X_Psttcr == 0)

dim(a)

# Randomly selecting pseudo-absences

n_l <- nrow(a)

n <- min(134, n_l)

l <- a[sample(1:n_l, n), ]

st_write(l, "grid_a.shp")

mpa1<-dplyr::select(l, X_EcP.18, X_Psttcr,X_Qldcht,X_Qlgrnd,X_Qlmltf,X_Qlprvf,X_Slvrtc,X_Schsrr ,X_Trsnvv,X_Thrpss,X_Vchysc,X_Vchyse,X_Vchysr,X_Vchysth,X_Vchystc)

mpa2<-dplyr::select(p, X_EcP.18, X_Psttcr,X_Qldcht,X_Qlgrnd,X_Qlmltf,X_Qlprvf,X_Slvrtc,X_Schsrr ,X_Trsnvv,X_Thrpss,X_Vchysc,X_Vchyse,X_Vchysr,X_Vchysth,X_Vchystc)

head(mpa1)
head(mpa2)

inc<-rbind(mpa1,mpa2)

head(inc)

dados_inc <- st_drop_geometry(inc)

head(dados_inc)

write.xlsx(dados_inc, "mpa.xlsx")


##THE END##

rm(list=ls())

