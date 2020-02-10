#-------------------------------------------------------------------------------------------------------------------
#this code makes figure 2
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
# First load and post process the output of STEPPS
library(rstan)
library(RColorBrewer)
library(fields)
library(sp)
library(maptools)
#-----------------------------------------------------------------------------------------------
repo_location <- wd
setwd(paste(repo_location,'prediction/',sep=''))
help.fun.loc <- 'utils/'
data.loc <- 'output_nb/'
plot.loc <- 'plots/'
#-----------------------------------------------------------------------------------------------
us.shp <- readShapeLines('data/map_data/us_alb.shp',proj4string=CRS('+init=epsg:3175'))
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Oak','Hemlock'))

K <- length(taxa)



######################################################################################################################
#load pre-processed output
######################################################################################################################
proc.output.loc <- 'process_output/'

mean.proportion <- readRDS(paste(proc.output.loc,'combined_mean.RDS',sep=''))  
median.proportion <-readRDS(paste(proc.output.loc,'combined_median.RDS',sep=''))
lb.proportion <- readRDS(paste(proc.output.loc,'combined_lb.RDS',sep=''))
ub.proportion <- readRDS(paste(proc.output.loc,'combined_ub.RDS',sep=''))
diff_bounds_proportion <- ub.proportion - lb.proportion


#----------------------------------------------------------------------------------------------------------------------
#load vegetation data and coordinates
#
#----------------------------------------------------------------------------------------------------------------------
#load vegetation data (contains median of vegetation)

load(paste0(repo_location,'calibration/data/for_prediction_110_sites_only_abies_new_species.RData'))
veg <- as.matrix(r)
#r = meidan of vegetation
# veg_coords = coordinates of vegetation 

#Transform pollen coordinates
sputm <- SpatialPoints(pollen_coords, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
pollen_coord_us <- spgeo@coords
#---------------------------------------------------------------------------------------------------------------------
#load REVEALS 
reveals_reconstruction <- readRDS(paste0(wd,'Reveals_NEUS/output/unc_stand_veg_pred_ppe_literature_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral_1e+05_regional_vegetation_aggregated.RDS'))
reveals_coord <- readRDS(paste0(wd,'Reveals_NEUS/output/unc_stand_ratio_pred_ppe_literature_fallspeed_constant_composition_larch_lake_size_correct_dwm_gpm neutral_1e+05_regional_vegetation_aggregated_coord.RDS'))

reveals_median <- reveals_reconstruction[,grep('median',colnames(reveals_reconstruction))] 
reveals_mean <- reveals_reconstruction[,grep('mean',colnames(reveals_reconstruction))] 
colnames(reveals_median)[grep('Other',colnames(reveals_median))] <- c('Other conifer','Other hardwood') 
colnames(reveals_mean)[grep('Other',colnames(reveals_mean))] <- c('Other conifer','Other hardwood') 

#Transform reveals coordinates in to US coordinates
sputm <- SpatialPoints(reveals_coord[,c('x','y')], proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
reveals_coord_us <- spgeo@coords


#---------------------------------------------------------------------------------------------------------------------
#find grid cells where REveals has a prediction
bounding_coords <- matrix(ncol=4, nrow = nrow(reveals_coord))
bounding_coords[,1] <- reveals_coord[,'x']-0.5
bounding_coords[,2] <- reveals_coord[,'x']+0.5
bounding_coords[,3] <- reveals_coord[,'y']-0.5
bounding_coords[,4] <- reveals_coord[,'y']+0.5

colnames(bounding_coords) <- c('lon.min','lon.max','lat.min','lat.max')

spgeo <- SpatialPoints(veg_coords, proj4string=CRS("+init=epsg:3175"))
sputm <- spTransform(spgeo, CRS("+init=epsg:4326"))
veg_coords_utm  <- sputm@coords #transformation is successful
colnames(veg_coords_utm) <- c('lon','lat')

coords_assign <-  matrix(ncol = 1,nrow =nrow(veg_coords))

for(x in 1:nrow(bounding_coords)){
  coords_assign[((veg_coords_utm[,'lon'] > bounding_coords[x,"lon.min"]) & (veg_coords_utm[,'lon'] < bounding_coords[x,"lon.max"]) &
                   (veg_coords_utm[,'lat'] > bounding_coords[x,"lat.min"]) & (veg_coords_utm[,'lat'] < bounding_coords[x,"lat.max"]))] <- x
}


#--------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------
#Figure begins
#----------------------------------------------------------------------------------------------------------------------
setwd(paste(repo_location,'/figures',sep=''))
plot.loc <- 'figures/'







taxa_plot <- c('Beech','Hemlock','Oak','Spruce')
colours <- rev(brewer.pal(10,'RdYlBu'))
east <- sort(unique(veg_coords[,'meters.east']))
north <- sort(unique(veg_coords[,'meters.north']))
breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)


#pdf(paste(plot.loc,'figure_2_all_dec.pdf',sep=''),height = 44,width = 10)
pdf(paste(plot.loc,'figure_2.pdf',sep=''),height = 16,width = 10)
#par(mfrow=c(11,3))
par(mfrow=c(4,3))
par(oma=c(1,2,2,2),mar=c(0.1,0.1,0.1,0.1))
for(y in taxa_plot){
par(cex=1)  

# plot vegetatin data  
  categories <- cut(veg[,y],breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]

  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  points(veg_coords,col=colours.plot,pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  mtext(side = 2,font=2,line=0,text = y,cex = 2)
  if(y==taxa_plot[1]) {
    mtext(side = 3,font=2,line=0,text = 'Township survey',cex = 2)
  }
  #-
#-------------------------------------------------------------------------------------------------------------------  
# plot STEPPS reconstruction 
  #categories <- cut(median.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  categories <- cut(mean.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]
  
  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  if(y%in%taxa_plot[seq(1,11,4)]){
    image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,add = TRUE,horizontal=TRUE,
             legend.width = 6,legend.mar = 2.5,legend.lab='Proportion')
  }
  points(veg_coords,col=colours.plot,pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  if(y==taxa_plot[1]) {
    mtext(side = 3,font=2,line=0,text = 'STEPPS',cex = 2)
  }
#-----------------------------------------------------------------------------------------------------------------  
# plot REVEALS reconstruction   
  #categories <- cut(reveals_median[,grep(y,colnames(reveals_median))]/100,breaks,include.lowest = TRUE,labels = FALSE)
  categories <- cut(reveals_mean[,grep(y,colnames(reveals_mean))]/100,breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]
  colours.plot.final <- matrix(ncol=1, nrow=length(coords_assign))
  
  for(i in 1:length(colours.plot)){
    colours.plot.final[coords_assign==i] <- colours.plot[i]
  }  
  colours.plot.final[is.na(colours.plot.final)] <- 'gray'
  
  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  points(veg_coords,col=colours.plot.final,pch = 15,cex=0.5)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  if(y%in%taxa_plot[seq(1,11,4)]) {
    mtext(side = 3,font=2,line=0,text = 'REVEALS',cex = 2)
  }
#-----------------------------------------------------------------------------------------------------------------  
}
dev.off()



###########################################################################################################################
#
###########################################################################################################################
taxa_plot <- taxa[!taxa%in%c('Chestnut','Tamarack')]
pdf(paste(plot.loc,'SF3.pdf',sep=''),height = 52,width = 13.3333)
par(mfrow=c(13,4))
par(oma=c(2,2,2,2),mar=c(0.1,0.1,0.1,0.1))
for(y in taxa){
  par(cex=1)  
  
  # plot vegetatin data  
  
  categories <- cut(median.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]
  
  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  points(veg_coords,col=colours.plot,pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  mtext(side = 2,font=2,line=0,text = y,cex = 2)
  if(y==taxa_plot[1]) {
    mtext(side = 3,font=2,line=0,text = 'Median',cex = 2)
  }
  categories <- cut(diff_bounds_proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]
  
  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  points(veg_coords,col=colours.plot,pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  if(y==taxa_plot[1]) {
    mtext(side = 3,font=2,line=0,text = 'Prediction Interval',cex = 2)
  }
  #-
  #-------------------------------------------------------------------------------------------------------------------  
  # plot STEPPS reconstruction 
  #categories <- cut(median.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  categories <- cut(lb.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]
  
  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,add = TRUE,horizontal=TRUE,
             legend.width = 4,legend.mar = 3.5)
  points(veg_coords,col=colours.plot,pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  if(y==taxa_plot[1]) {
    mtext(side = 3,font=2,line=0,text = 'Lower Bound',cex = 2)
  }
  #-----------------------------------------------------------------------------------------------------------------  
  # plot REVEALS reconstruction   
  #categories <- cut(reveals_median[,grep(y,colnames(reveals_median))]/100,breaks,include.lowest = TRUE,labels = FALSE)
  categories <- cut(ub.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]
  colours.plot.final <- matrix(ncol=1, nrow=length(coords_assign))
  
  for(i in 1:length(colours.plot)){
    colours.plot.final[coords_assign==i] <- colours.plot[i]
  }  
  colours.plot.final[is.na(colours.plot.final)] <- 'gray'
  
  plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north))  
  points(veg_coords,col=colours.plot,pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 1)
  plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
  if(y==taxa_plot[1]) {
    mtext(side = 3,font=2,line=0,text = 'Upper Bound',cex = 2)
  }
  #-----------------------------------------------------------------------------------------------------------------  
  # plot REVEALS reconstruction with stepps ppe   
  # categories <- cut(reveals_median_stepps_ppe[,grep(y,colnames(reveals_median_stepps_ppe))]/100,breaks,include.lowest = TRUE,labels = FALSE)
  # colours.plot <- colours[categories]
  # 
  # plot(us.shp[c(3,9,12,14,16:20),],xlim=range(east),ylim=range(north),main =ifelse(y==taxa_plot[1],'REVEALS STEPPS PPE',''),cex.main = 2)  
  # points(reveals_coord_us,col=colours.plot,pch = 15,cex =4)
  # points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
  # plot(us.shp[c(3,9,12,14,16:20),],add=TRUE)
}
dev.off()

