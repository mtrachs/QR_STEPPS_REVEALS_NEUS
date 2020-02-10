##########################################################################################################################
# This script is used to compare vegetation reconstructions and vegetation observations using the chord distance 
##########################################################################################################################


#load data
library(rstan)
library(RColorBrewer)
library(fields)
library(sp)
library(maptools)
library(colorRamps)
library(rioja)
#-----------------------------------------------------------------------------------------------
setwd(paste0(wd,'prediction/'))
help.fun.loc <- 'utils/'
data.loc <- 'output_nb//'
plot.loc <- 'plots/'
#-----------------------------------------------------------------------------------------------
us.shp <- readShapeLines('data/map_data/us_alb.shp',proj4string=CRS('+init=epsg:3175'))
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Oak','Hemlock'))
K <- length(taxa)

proc.output.loc <- 'process_output/'

mean.proportion <- readRDS(paste(proc.output.loc,'combined_mean.RDS',sep=''))  
median.proportion <-readRDS(paste(proc.output.loc,'combined_median.RDS',sep=''))
lb.proportion <- readRDS(paste(proc.output.loc,'combined_lb.RDS',sep=''))
ub.proportion <- readRDS(paste(proc.output.loc,'combined_ub.RDS',sep=''))

#--------------------------------------------------------------------------------------------------------------------------
#load vegetation data
load(paste0(wd,'calibration/data/for_prediction_110_sites_only_abies_new_species.RData'))
r <- readRDS(paste0(wd,'calibration/data/veg_mean_dissimilarity.RDS'))

########################################################################################################################

mean.rec.sqrt <- sqrt(mean.proportion)
mean.obs.sqrt <- sqrt(r)

#-------------------------------------------------------------------------------
core.loc.rec <- mean.rec.sqrt[unique(idx_cores),]
core.loc.obs <- mean.obs.sqrt[unique(idx_cores),]

dissimilarity.dep.sites <- 
sapply(1:nrow(core.loc.rec),function(x){
  as.matrix(paldist((rbind(core.loc.rec[x,],core.loc.obs[x,])),dist.method='euclidean'))[1,2]
})

#--------------------------------------------------------------------------------
other.loc.rec <- mean.rec.sqrt[!(1:nrow(mean.rec.sqrt))%in%unique(idx_cores),]
other.loc.obs <- mean.obs.sqrt[!(1:nrow(mean.obs.sqrt))%in%unique(idx_cores),]

dissimilarity.other.sites <- 
  sapply(1:nrow(other.loc.rec),function(x){
    as.matrix(paldist((rbind(other.loc.rec[x,],other.loc.obs[x,])),dist.method='euclidean'))[1,2]
  })


#----------------------------------------------------------------------------------
#compare means
mean.other <- mean(dissimilarity.other.sites)
mean.dep <- mean(dissimilarity.dep.sites)

t.test(dissimilarity.other.sites^2,dissimilarity.dep.sites^2)

####################################################################################################################
#dissimilarites of all sites
####################################################################################################################
# dissimilarity.tot <- 
#   sapply(1:nrow(mean.rec.sqrt),function(x){
#     as.matrix(paldist((rbind(mean.rec.sqrt[x,],mean.obs.sqrt[x,])),dist.method='euclidean'))[1,2]
#   })
# 
# dissimilarity.tot <- dissimilarity.tot^2


dissimilarity.tot <- 
  sapply(1:nrow(mean.rec.sqrt),function(x){
    as.matrix(paldist(rbind(mean.proportion[x,],r[x,])))[1,2]
  })

#####################################################################################################################
# plot dissimilarities
#####################################################################################################################
breaks <- 2*c(0,0.025,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
categories <- cut(dissimilarity.tot,breaks,include.lowest = TRUE,labels = FALSE)
colours <- rev(brewer.pal(10,'RdYlBu'))

colours.plot <- colours[categories]

east <- sort(unique(veg_coords[,'meters.east']))
north <- sort(unique(veg_coords[,'meters.north']))

breaks1 <- 2*c(0,0.025,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
breaks2 <- 2*c(0,0.025,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)

sputm <- SpatialPoints(pollen_coords, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
pollen_coord_us <- spgeo@coords

x1 <- '110_sites_umw'#strsplit(x,'.csv')[[1]][1]

pdf(paste(wd,'/figures/figures/Figure4.pdf',sep=''),width = 6,height = 5)
  par(oma=c(1,1,1,2),cex = 0.8,cex.axis = 0.8)

  plot(us.shp,xlim=range(east),ylim=range(north),main= "Dissimilarity Map")
  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
           breaks= breaks1,lab.breaks = breaks2,cex.axis = 0.8,add=TRUE,legend.lab = 'Sq. Chord Distance',legend.line = 2.5)
  points(veg_coords,col=colours.plot,pch = 15)
  plot(us.shp,add=TRUE)
  points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
dev.off()
#})

#####################################################################################################################
# check for analogue quality
#####################################################################################################################
# library(analogue)
# 
# data("Pollen")
# Pollen[is.na(Pollen)] <- 0
# pollen.perc <- Pollen/rowSums(Pollen)
# sq.chord.dist.pollen <- paldist(pollen.perc)
# # Richard Telford's rule of Thumb (might be problematic with large calibration datasets)
# quantile(sq.chord.dist.pollen,probs=c(0.05,0.1))
# 
# #mcarlo.dist <- mcarlo(pollen.perc,nsamp=100,typ='bootstrap',method='SQchord')
# 
# 
# library(palaeoSig)
# data("arctic.pollen")
# 
# dist.arctic.pollen <- paldist(arctic.pollen/100)
# # Richard Telford's rule of Thumb (might be problematic with large calibration datasets)
# quantile(dist.arctic.pollen,probs=c(0.05,0.1))



#---------------------------------------------------------------------------------------------------------------------
######################################################################################################################
#estimate dissimilarity between REVEALS and Paciorek et al. 
######################################################################################################################

#---------------------------------------------------------------------------------------------------------------------
#load REVEALS predictions
reveals_gridded <- readRDS(paste0(wd,'Reveals_NEUS/data/reveals_mean_gridded.RDS'))

colnames(reveals_gridded) <- unlist(strsplit(colnames(reveals_gridded),'[.]'))[seq(1,2*ncol(reveals_gridded),2)]

colnames(reveals_gridded)[grep('Other',colnames(reveals_gridded))] <- c('Other conifer','Other hardwood') 
#again load meanof Paciorek et al. 
r <- readRDS(paste0(wd,'calibration/data/veg_mean_dissimilarity.RDS'))

r_reveals <- r[,colnames(r)%in%colnames(reveals_gridded)]
r_reveals <- r_reveals[!is.na(reveals_gridded[,"Beech"]),]
r_reveals <- r_reveals/rowSums(r_reveals)
reveals_gridded_complete <- reveals_gridded[!is.na(reveals_gridded[,"Beech"]),] 

dissimilarity.reveals <- 
  sapply(1:nrow(r_reveals),function(x){
    as.matrix(paldist(rbind(r_reveals[x,],reveals_gridded_complete[x,])))[1,2]
  })


#---------------------------------------------------------------------------------------------------------------------
stepps_11_taxa <- mean.proportion[,colnames(mean.proportion)%in%colnames(reveals_gridded)]
stepps_11_taxa <- stepps_11_taxa/rowSums(stepps_11_taxa)
stepps_11_taxa_coverage_reveals <- stepps_11_taxa[!is.na(reveals_gridded[,"Beech"]),]
r_11_taxa <- r[,colnames(r)%in%colnames(reveals_gridded)]
r_11_taxa <- r_11_taxa/rowSums(r_11_taxa)

dissimilarity.stepps_11_taxa<- 
  sapply(1:nrow(r_11_taxa),function(x){
    as.matrix(paldist(rbind(r_11_taxa [x,],stepps_11_taxa[x,])))[1,2]
  })

t.test(dissimilarity.reveals,dissimilarity.stepps_11_taxa)

dissimilarity.stepps_11_taxa_coverage_reveals<- 
  sapply(1:nrow(r_reveals),function(x){
    as.matrix(paldist(rbind(r_reveals[x,],stepps_11_taxa_coverage_reveals[x,])))[1,2]
  })

t.test(dissimilarity.reveals,dissimilarity.stepps_11_taxa_coverage_reveals,paired=TRUE)



######################################################################################################################
#
######################################################################################################################
comparison_dissimilarities <- cbind(summary(dissimilarity.tot),
      summary(dissimilarity.stepps_11_taxa_coverage_reveals),
      summary(dissimilarity.reveals))

write.csv(round(comparison_dissimilarities,3),paste(data.loc,'comparison_dissimilarities.csv',sep=''))


