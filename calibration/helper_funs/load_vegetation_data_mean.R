library(ncdf4)
library(stepps)
veg <- nc_open(paste0(wd,'calibration/data/SetTreeComp_Level2_v1.0.nc'))
veg.names <- names(veg$var) #gives the names of the genera available 
#load a taxon at a time
veg.comp <- lapply(veg.names,function(x) {
  ncvar_get(veg,varid=x)
})

names(veg.comp) <- veg.names

#find coordinates
meters.east <- veg$dim$x$vals
meters.north <- veg$dim$y$vals

# sample is a sample of a posterior distribution of 250 samples that is retained

#plot a pdf of a posterior
plot(density(veg.comp$Oak[50,50,]))


#----------------------------------------------------------------------------------------------------------------------
#prepare that data fro STEPPS
#Simon uses structure x coordinate y coordinate,region, water, and then vegetation

#first let's prepare
#----------------------------------------------------------------------------------------------------------------------
veg.mean.post <- 
  sapply(names(veg.comp),function(z) {
    apply(veg.comp[[z]],c(1,2),function(x) mean(x,na.rm=TRUE))
  })

total.dat <- as.data.frame(cbind(rep(meters.east,length(meters.north)),rep(meters.north,each = length(meters.east)),
                                 veg.mean.post))
colnames(total.dat)[1:2] <-c('meters.east','meters.north')

#source(paste(help.fun.loc,'load_vegetation_data_median.R',sep=''))
veg_mean <- total.dat
colnames(veg_mean)[1:2] <- c('meters.east','meters.north')

#veg_mean <- readr::read_csv('data/composition_v0.3.csv')
coords.neus <- matrix(ncol=2,c(rep(c(-79,-67),each =2),rep(c(39.5,47.5),2))) 
coords.neus <- as.data.frame(coords.neus)
colnames(coords.neus) <- c('Lon','Lat')

pol_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                     '+init=epsg:4326', 
                     '+init=epsg:4326')

veg_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                     '+init=epsg:4326', 
                     '+init=epsg:3175')

#extract vegetation data that is in the bounding box
veg_mean <- veg_mean[(veg_mean$meters.east>veg_box[1]) & (veg_mean$meters.north>veg_box[2]),]



veg_mean <- veg_mean[!is.na(veg_mean$Ash),]
veg_mean <- replace(veg_mean,veg_mean==0,2e-7) #does this have an effect

setwd(paste0(wd,'calibration/'))
data.loc <- 'data/'

veg_table <- readr::read_csv(paste(data.loc,'/veg_trans_edited_only_abies.csv',sep=''))
veg_trans <- translate_taxa(veg_mean, veg_table ,id_cols = colnames(veg_mean)[1:2])

r <- veg_trans[,!colnames(veg_trans)%in%c("meters.east","meters.north")]
r[,"Other hardwood"] <- rowSums(r[,c("Other hardwood","Elm","Poplar")])
r <- r[,!colnames(r)%in%c("Elm","Poplar")]

saveRDS(r,paste(data.loc,'veg_mean_dissimilarity.RDS',sep=''))
