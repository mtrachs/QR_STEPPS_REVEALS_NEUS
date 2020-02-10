library(raster)
library(sp)
library(disqover)
library(reshape2)
library(ggplot2)
library(dplyr)
library(neotoma)
library(maps)
library(RColorBrewer)
library(fields)
library(zipfR)
library(rioja)
###################################################################################################################################
#load pollen data
 
setwd(paste0(wd,'Reveals_NEUS'))
data.loc <-'data/'
plot.loc <-'figures/'

load(paste0(wd,'calibration/data/elicitation_neus_certainty_median_110_sites_only_abies_new_species.RData'))

y1 <- y[,!colnames(y)%in%'Chestnut']
#reveals wants and age estimate in column1 (or at least no pollen data)
y2 <- cbind(1:nrow(y1),y1)


taxa <- colnames(y1)

###################################################################################################################################
####################################################################################################
#read REVEALS estimates
####################################################################################################
svs = read.csv(paste(data.loc,'svs_LC6K.csv',sep=''), sep=',', header=TRUE, stringsAsFactors=FALSE)

svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)
svs_agg$taxon[svs_agg$taxon=='Alder'] <- 'Other hardwood'
svs_agg$taxon[svs_agg$taxon=='Fir'] <- 'Other conifer'
svs_agg$taxon[svs_agg$taxon=='Larch'] <- 'Tamarack'
svs_agg <- svs_agg[svs_agg$taxon%in%taxa,]
svs_agg <- svs_agg[order(svs_agg$taxon),]

ppes <- readRDS(paste(data.loc,'PPEs_agg.RDS',sep=''))
ppes$taxon <- as.character(ppes$taxon)
ppes$taxon[ppes$taxon=='Alder'] <- 'Other hardwood'
ppes$taxon[ppes$taxon=='Fir'] <- 'Other conifer'
ppes$taxon[ppes$taxon=='Larch'] <- 'Tamarack'
ppes <- ppes[ppes$taxon%in%taxa,]
ppes <- ppes[order(ppes$taxon),]

params <- cbind(ppes[,c(1,3:4)],svs_agg$sv)
colnames(params) <- c('taxon','ppe','ppe.error','fallspeed') 
#################################################################################################
#test effects of lake radii
rad.lake <- c(5.5,55,550,5500)
recon.lake.size <-
lapply(rad.lake,function(radius){
a <- REVEALSinR(pollen = y2,#"~/Reveals_NEUS/data/reveals_test.csv",
                params = params,#"~/Reveals_NEUS/data/reveals_input_params.csv",
                dwm        = 'GPM neutral',
                tBasin     = "lake",
                dBasin     = 2*radius, # diameter!
                regionCutoff = 100000,#radius,
                n      = 1000)

recon.reveals <- a[,grep('mean',colnames(a))]
})

#numbers reported in paper
effect_lake_size  <- sapply(1:nrow(recon.lake.size[[1]]),function(x){
                              paldist2(recon.lake.size[[1]][x,]/100,recon.lake.size[[4]][x,]/100)
})

summary(effect_lake_size)

################################################################################################################################
#test effects of regional vegetation
rad.veg <- c(50000,100000,200000,300000,400000)
recon.tot <-
  lapply(rad.veg,function(radius){
    a <- REVEALSinR(pollen = y2,
                    params = params,
                    dwm        = 'GPM neutral',
                    tBasin     = "lake",
                    dBasin     = 2*550, # diameter!
                    regionCutoff = radius,
                    n      = 1000)
    
    recon.reveals <- a[,grep('mean',colnames(a))]
  })


effect_veg_rad  <- sapply(1:nrow(recon.tot[[1]]),function(x){
  paldist2(recon.tot[[1]][x,]/100,recon.tot[[4]][x,]/100)
})

summary(effect_veg_rad)

#################################################################################################################
# extract K from Martin Theurkauf's code
source(paste0('R/REVEALS_K.R'))
files.reveals <- list.files(paste0(wd,'disqover/R'))
sapply(files.reveals[grep('.R',files.reveals)],function(x) source(paste0(wd,'disqover/R/',x)))

rad.lake <- c(5.5,55,550,5500)
K <-
  sapply(rad.lake,function(radius){
    a <- REVEALSinR_K(pollen = y2[1,],
                    params = params,
                    dwm        = 'GPM neutral',
                    tBasin     = "lake",
                    dBasin     = 2*radius, # diameter!
                    regionCutoff = 100000,#radius,
                    n      = 1000)
    
    a[[2]]
  })

rownames(K) <- ppes$taxon
colnames(K) <- rad.lake

saveRDS(K,paste(data.loc,'K.RDS',sep=''))

