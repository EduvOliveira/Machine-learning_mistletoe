#' ---
# Disentangling biotic and abiotic drivers of a neotropical mistletoe 
# By Eduardo V. S. Oliveira
# 28/01/2026
#' ---

## Niche equivalency and similarity tests

install.packages(c('raster','virtualspecies','ecospat', 'ade4'))

library(raster)
library(virtualspecies)
library(ecospat)
library(ade4)

# Load bioclimatic data

path<-setwd("J:/Eduardo_Vinícius/R/ParasitaDispersorModelagemNicho/Reanalise/sources")

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
bio4<-preds_AS$bio4
bio5<-preds_AS$bio5
bio12<-preds_AS$bio12
bio14<-preds_AS$bio14
bio15<-preds_AS$bio15
bio18<-preds_AS$bio18
bio19<-preds_AS$bio19

my_preds<-stack(bio2,bio4,bio5,bio12,bio14,bio15,bio18,bio19)
plot(my_preds)

# Entering data on species occurrences

occ<-read.csv("occ_clean_rea.csv", h=TRUE)

occ<-occ[,2:4]

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

# Selection of plant species

levels(as.factor(occ$spp))

sp.choice<- c('e_cristata','p_robustus','q_dichotoma', 'q_grandiflora','q_multiflora','q_parviflora','s_convallariodora','s_ruficapillus','t_sayaca','t_viridis','v_cinnamomea','v_elliptica','v_rufa','v_thyrsoidea','v_tucanorum') 

sp.combn<-combn(sp.choice,2) 
nsp<-15


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
