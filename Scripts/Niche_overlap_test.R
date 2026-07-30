#' ---
# Disentangling biotic and abiotic drivers of a neotropical mistletoe 
# By Eduardo V. S. Oliveira
# 28/07/2026
#' ---

## Niche equivalency and similarity tests

install.packages(c('raster','virtualspecies','ecospat', 'ade4'))

library(raster)
library(virtualspecies)
library(ecospat)
library(ade4)

# Load bioclimatic data

path<-setwd("J:/path_to_the_folder")

lst <- list.files(path=path,pattern='tif$',full.names = T) 
preds<-stack(lst)
names(preds)
names(preds) <- paste0("bio", 1:19)

# Cutting to South America
Am_Sul<-shapefile(choose.files())  

A_Sul<-aggregate(Am_Sul,dissolve=TRUE)
plot(A_Sul) 

ame<-shapefile(choose.files()) 
ame

preds_AS<-list()
for(i in 1:19) preds_AS[[i]]<-mask(preds[[i]],ame)

preds_AS<-stack(preds_AS)
preds_AS<-crop(preds_AS, c(-85, -30, -60, 15))

plot(preds_AS[[1]])
plot(ame, add=T)

# checking collinearity

var<-removeCollinearity(preds_AS, multicollinearity.cutoff = 0.7, select.variables = TRUE, sample.points = FALSE, plot = FALSE)

# Selecting the variables

bio2<-preds_AS$bio2
bio3<-preds_AS$bio3
bio10<-preds_AS$bio10
bio13<-preds_AS$bio13
bio15<-preds_AS$bio15
bio18<-preds_AS$bio18
bio19<-preds_AS$bio19

my_preds<-stack(bio2,bio4,bio5,bio12,bio14,bio15,bio18,bio19)
plot(my_preds)

# Entering data on species occurrences

occ<-read.csv("occ_clean.csv", h=TRUE)

occ<-occ[,1:3]

occ$species <- gsub(" ", "_", occ$species)

# Removes species with fewer than 5 records

which(table(occ$species) < 5)

sp.validas <- names(which(table(occ$species) >= 5))
occ <- occ[occ$species %in% sp.validas, ]

# Extracting climate data from Raster
clim.bkg<-na.exclude(data.frame(extract(my_preds,ame)))
clim.occ<-na.exclude(data.frame(extract(my_preds,occ[,2:3])))

# Species list
sp.list<-levels(as.factor(occ[,1]))

sp.nbocc<-c()

for (i in 1:length(sp.list)){
  sp.nbocc<-c(sp.nbocc,length(which(occ[,1] == sp.list[i])))
  
} #Calculate the number of observed occurrences
nb.sp<-length(sp.list) 

# PCA-env calibration
pca.env <-dudi.pca(clim.bkg, center = T, scale = T, scannf = F, nf = 2)

# Selection of species

levels(as.factor(occ$spp))

sp.choice<- c('Campomanesia_adamantium','Diplusodon_hirsutus', 'Elaenia_cristata','Miconia_albicans','Miconia_ferruginata','Psittacanthus_robustus','Qualea_cordata','Qualea_dichotoma','Qualea_grandiflora','Qualea_multiflora','Qualea_parviflora','Salvertia_convallariodora','Schistochlamys_ruficapillus','Tersina_viridis','Thraupis_sayaca','Trembleya_laniflora','Vochysia_cinnamomea','Vochysia_elliptica','Vochysia_rufa','Vochysia_thyrsoidea','Vochysia_tucanorum')  

sp.combn<-combn(sp.choice,2) 
nsp<-22

# Arrays for storing
overlap<-matrix(nrow=nsp,ncol=nsp,dimnames = list(sp.choice,sp.choice))		#Records overlap values

equivalency<-matrix(nrow=nsp,ncol=nsp,dimnames = list(sp.choice,sp.choice))	#records p-values for equivalence test

similarity<-matrix(nrow=nsp,ncol=nsp,dimnames = list(sp.choice,sp.choice))	#records p-values for similarity test

# Loop of niche quantification for each combination of species

for(i in 1:ncol(sp.combn)) {  
  spa<-sp.combn[1,i] 
  spb<-sp.combn[2,i] 
  clim.spa<-clim.occ[occ$spp==spa,] #env data for species a
  clim.spb<-clim.occ[occ$spp==spb,] #env data for species b
  
  # PCA scores
  scores.bkg<- pca.env$li	#scores for global climate
  scores.spa<- suprow(pca.env,clim.spa)$lisup				#scores for spa
  scores.spb<- suprow(pca.env,clim.spb)$lisup				#scores for spb
  
  # calculation of occurence density
  za<- ecospat.grid.clim.dyn(scores.bkg,scores.bkg,scores.spa,100)
  zb<- ecospat.grid.clim.dyn(scores.bkg,scores.bkg,scores.spb,100)
  
  # overlap corrected by availabilty of background conditions
  ecospat.niche.overlap(za,zb,cor=F) 
  
  # test of niche equivalency and similarity
  equ<-ecospat.niche.equivalency.test(za,zb,rep=100) 
  
  sim<-ecospat.niche.similarity.test(za,zb,rep=100,rand.type = 1) 
  
  
  #storage of values
  overlap[sp.combn[1,i],sp.combn[2,i]]<-ecospat.niche.overlap(za,zb,cor=T)[[1]] 	#store overlap value
  equivalency[sp.combn[1,i],sp.combn[2,i]]<-equ$p.D					#store equivalency value
  similarity[sp.combn[1,i],sp.combn[2,i]]<-sim$p.D				#store similarity 21 value
  
  #plot			
  ecospat.plot.niche(za,title=spa,name.axis1="PC1",name.axis2="PC2") 
  ecospat.plot.niche(zb,title=spb,name.axis1="PC1",name.axis2="PC2")
  ecospat.plot.contrib(pca.env$co,pca.env$eig)
  ecospat.plot.overlap.test(equ,"D","Equivalency")
  ecospat.plot.overlap.test(sim,"D","Similarity")
  
  #counter
  cat(paste(i))
}
########################################################
# This script was based on the codes provided in the paper Broennimann et al. 2012. Measuring ecological overlap from occurrence and spatial environmental data. Global Ecology and Biogeography 21:481-497.
########################################################

##THE END##

rm(list=ls())
