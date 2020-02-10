#-----------------------------------------------------------------------------------------------------------
#Figure 3 of manuscript
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
library(rstan)
library(RColorBrewer)
library(fields)
library(sp)
library(maptools)
library(rioja)
#------------------------------------------------------------------------------------------------------------------------
#compare REVEALS and STEPPS
#------------------------------------------------------------------------------------------------------------------------
reveals.loc <- paste0(wd,'Reveals_NEUS/')
reveals.data <- 'output/'
help.fun.loc <- paste0(wd,'figures/utils/')
source(paste(help.fun.loc,'pairs_functions.R',sep=''))


reveals.output <- readRDS(paste(reveals.loc,reveals.data,
                                'unc_stand_veg_pred_ppe_literature_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral_1e+05_regional_vegetation.RDS',sep=''))
#'veg_pred_ppe_stepps_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral.RDS',sep=''))

mean.sim.index <- grep('mean',colnames(reveals.output))
reveals.reconstructions <- reveals.output[,mean.sim.index]
colnames(reveals.reconstructions)[grep('Other',colnames(reveals.reconstructions))] <- c('Other conifer',"Other hardwood")
#is currently necessary because REVEALS stil simulates Elm and Poplar
reveals.reconstructions <- reveals.reconstructions[,!colnames(reveals.reconstructions)%in%c('Elm.meansim','Poplar.meansim')]
#reveals.coords <- reveals.output[,c()]
  
# sputm <- SpatialPoints(reveals.coords, proj4string=CRS("+init=epsg:4326"))
# spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
# reveals.coords.us <- spgeo@coords

#------------------------------------------------------------------------------------------------------------------------
#load vegetation data
#------------------------------------------------------------------------------------------------------------------------
load(paste0(wd,'calibration/data/for_prediction_110_sites_only_abies_new_species.RData'))
r <- readRDS(paste0(wd,'calibration/data/veg_mean_dissimilarity.RDS'))
veg <- as.matrix(r)
#veg_coords

coords.stepps <- pollen_coords
#-----------------------------------------------------------------------------------------------------------------------
#reduce reveals to 73 sites
#-----------------------------------------------------------------------------------------------------------------------
reveals.coords.us <- coords.stepps


#----------------------------------------------------------------------------------------------------------------------
#load vegetation predictions
#-----------------------------------------------------------------------------------------------
setwd(paste0(wd,'prediction/'))
help.fun.loc <- 'utils/'
data.loc <- 'output_nb/'
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





rec.stepps <- mean.proportion[,!(colnames(mean.proportion)%in%c('Chestnut','Tamarack'))]

#-----------------------------------------------------------------------------------------------------------------------
veg_dep_sites <- as.data.frame(veg[idx_cores,!(colnames(veg)%in%c('Chestnut','Tamarack','Poplar','Elm'))])
stepps_dep_sites <- as.data.frame(rec.stepps[idx_cores,])
reveals.reconstructions <- as.data.frame(reveals.reconstructions/100)
colnames(reveals.reconstructions) = colnames(veg_dep_sites)


#-------------------------------------------------------------------------------------------------------------------
#test correlations between REVEALS and STEPPS

linear_relation <- 
sapply(colnames(stepps_dep_sites),function(x){
  cors <- cor.test(stepps_dep_sites[,x],reveals.reconstructions[,x])
  data.frame(r = round(cors$estimate,2),p = round(cors$p.value,4))
})

monotonic_relation <- 
  sapply(colnames(stepps_dep_sites),function(x){
    cors <- cor.test(stepps_dep_sites[,x],reveals.reconstructions[,x],method='s')
    data.frame(r = round(cors$estimate,2),p = round(cors$p.value,4))
  })

#-------------------------------------------------------------------------------------------------------------
#plotselected taxa
main.taxa <- c('Beech','Hemlock','Oak','Spruce')
all_taxa_plot <- taxa[!taxa%in%c('Chestnut','Tamarack')]

#pdf('~/workflow_stepps_calibration/figures_manuscript/plots/Figure_3_all_taxa.pdf',height = 13.3333,width = 10)
pdf(paste0(wd,'figures/figures/Figure3.pdf'),height = 8,width = 7)
par(mfrow=c(4,3),oma=c(2,2,2,0.5),mar=c(1.5,1.5,0,0),cex.axis=1.25)
for(i in main.taxa) {
#for(i in all_taxa_plot) {  
     taxon <- data.frame(veg = veg_dep_sites[,i],STEPPS = stepps_dep_sites[,i],
                       REVEALS=reveals.reconstructions[,i])
  # plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),main = i,xlab='',ylab='')
  # points(taxon$STEPPS,taxon$veg,col=2)
  # points(taxon$REVEALS,taxon$veg,col=4)
  # points(taxon$STEPPS,taxon$REVEALS,col='orange')
  plot(taxon$STEPPS,taxon$veg,col=1,pch = 16,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  if(i==main.taxa[4]){
    mtext(side = 1,line=2.2,font = 2, text = 'STEPPS')
    axis(1)
  }
  mtext(side = 2,line=2.2,font = 2, text = 'TPS')
  mtext(side=3,line=0,text = i,font=2)
  abline(a=0,b=1,lty=2)
  axis(2)
  box()
  plot(taxon$REVEALS,taxon$veg,col=1,pch = 16,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  if(i==main.taxa[4]){
    mtext(side = 1,line=2.2,font = 2, text = 'REVEALS')
    axis(1)
  }
  mtext(side = 2,line=0,font = 2, text = 'TPS')
  mtext(side=3,line=0,text = i,font=2)
  abline(a=0,b=1,lty=2)
  box()
  plot(taxon$STEPPS,taxon$REVEALS,col=1,pch = 16,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  mtext(side = 2,line=0,font = 2, text = 'REVEALS')
  if(i==main.taxa[4]){
    mtext(side = 1,line=2.2,font = 2, text = 'STEPPS')
    axis(1)
  }
  mtext(side=3,line=0,text = i,font=2)
  abline(a=0,b=1,lty=2)
  box()
}
dev.off()


###########################################################################################################################
# similar plot for all taxa
############################################################################################################################


pdf(paste0(wd,'figures/figures/SF4.pdf'),height = 13.3333,width = 10)
par(mfrow=c(4,3),oma=c(2,2,2,0.5),mar=c(1.5,1.5,0,0),cex.axis=1.25)
for(i in all_taxa_plot) {  
  taxon <- data.frame(veg = veg_dep_sites[,i],STEPPS = stepps_dep_sites[,i],
                      REVEALS=reveals.reconstructions[,i])
  # plot(1,1,type='n',xlim=c(0,1),ylim=c(0,1),main = i,xlab='',ylab='')
  # points(taxon$STEPPS,taxon$veg,col=2)
  # points(taxon$REVEALS,taxon$veg,col=4)
  # points(taxon$STEPPS,taxon$REVEALS,col='orange')
  plot(taxon$STEPPS,taxon$veg,col=1,pch = 16,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  if(i%in%all_taxa_plot[c(4,8,11)]){
    mtext(side = 1,line=2.2,font = 2, text = 'STEPPS')
    axis(1)
  }
  mtext(side = 2,line=2.2,font = 2, text = 'TPS')
  mtext(side=3,line=0,text = i,font=2)
  abline(a=0,b=1,lty=2)
  axis(2)
  box()
  plot(taxon$REVEALS,taxon$veg,col=1,pch = 16,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  if(i%in%all_taxa_plot[c(4,8,11)]){
    mtext(side = 1,line=2.2,font = 2, text = 'REVEALS')
    axis(1)
  }
  mtext(side = 2,line=0,font = 2, text = 'TPS')
  mtext(side=3,line=0,text = i,font=2)
  abline(a=0,b=1,lty=2)
  box()
  plot(taxon$STEPPS,taxon$REVEALS,col=1,pch = 16,xlim=c(0,1),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  mtext(side = 2,line=0,font = 2, text = 'REVEALS')
  if(i%in%all_taxa_plot[c(4,8,11)]){
    mtext(side = 1,line=2.2,font = 2, text = 'STEPPS')
    axis(1)
  }
  mtext(side=3,line=0,text = i,font=2)
  abline(a=0,b=1,lty=2)
  box()
}
dev.off()

