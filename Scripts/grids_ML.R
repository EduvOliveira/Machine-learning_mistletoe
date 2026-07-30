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
library(exactextractr)

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

mpa<-dplyr::select(g, g, X_Cmpmn_, X_Dplsd_,X_Eln_cr,X_Mcn_lb,X_Mcn_fr,X_Psttc_,X_Ql_crd,X_Ql_dch,X_Ql_grn,X_Ql_mlt,X_Ql_prv,X_Slvrt_,X_Schst_,X_Trsn_v,X_Thrps_,X_Trmbl_,X_Vchys_c,X_Vchys_l,X_Vchys_r,X_Vchys_th,X_Vchys_tc)

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

mpa1<-dplyr::select(l, X_Cmpmn_, X_Dplsd_,X_Eln_cr,X_Mcn_lb,X_Mcn_fr,X_Psttc_,X_Ql_crd,X_Ql_dch,X_Ql_grn,X_Ql_mlt,X_Ql_prv,X_Slvrt_,X_Schst_,X_Trsn_v,X_Thrps_,X_Trmbl_,X_Vchys_c,X_Vchys_l,X_Vchys_r,X_Vchys_th,X_Vchys_tc)

mpa2<-dplyr::select(p, X_Cmpmn_, X_Dplsd_,X_Eln_cr,X_Mcn_lb,X_Mcn_fr,X_Psttc_,X_Ql_crd,X_Ql_dch,X_Ql_grn,X_Ql_mlt,X_Ql_prv,X_Slvrt_,X_Schst_,X_Trsn_v,X_Thrps_,X_Trmbl_,X_Vchys_c,X_Vchys_l,X_Vchys_r,X_Vchys_th,X_Vchys_tc)

head(mpa1)
head(mpa2)

inc<-rbind(mpa1,mpa2)

head(inc)

dados_inc <- st_drop_geometry(inc)

head(dados_inc)

write.xlsx(dados_inc, "mpa.xlsx")

# Inputting the variables

setwd("J:/path_to_the_folder")

lst <- list.files(path=".",pattern='tif$',full.names = T) 
preds<-stack(lst)
names(preds)

# Selecting the independent variables

bio2<-preds$wc2.1_5m_bio_02
bio3<-preds$wc2.1_5m_bio_03
bio10<-preds$wc2.1_5m_bio_10
bio13<-preds$wc2.1_5m_bio_13
bio15<-preds$wc2.1_5m_bio_15
bio18<-preds$wc2.1_5m_bio_18
bio19<-preds$wc2.1_5m_bio_19

my_preds<-stack(bio2,bio3,bio10,bio13,bio15,bio18,bio19)

plot(my_preds)

d1<-data.frame(exact_extract(preds, l, 'mean'))#plots

d2<-data.frame(exact_extract(preds, p, 'mean'))#plots

head(d1)
head(d2)

names(d1) <- gsub("mean.wc2.1_5m_", "", names(d1))
names(d2) <- gsub("mean.wc2.1_5m_", "", names(d2))

d<-rbind(d1,d2)

head(d)

# Organizing the presence-absence matrix

a<-shapefile(here("South_America.shp"))

mpa<-read.xlsx(here("resub","mpa.xlsx"))

dmpa<-cbind(mpa,d)

head(dmpa)

names(dmpa) <- gsub("X_", "", names(dmpa))

dmpa<-dmpa %>% relocate(Psttc_)

names(dmpa)

dmpa<-dmpa %>% rename(Probu= Psttc_, Cadam=Cmpmn_,Dhirs=Dplsd_,Ecris=Eln_cr,Malbi=Mcn_lb, Mferr=Mcn_fr, Qcord=Ql_crd,Qdich=Ql_dch,Qgran=Ql_grn,Qmult=Ql_mlt,Qparv=Ql_prv,Sconv=Slvrt_,Srufi=Schst_,Tviri=Trsn_v, Tsaya=Thrps_, Tlani=Trmbl_,Vcinn=Vchys_c,Velli=Vchys_l,Vrufa=Vchys_r,Vthyr=Vchys_th,Vtuca=Vchys_tc)

write.xlsx(dmpa, "mpa_data.xlsx")


##THE END##

rm(list=ls())

