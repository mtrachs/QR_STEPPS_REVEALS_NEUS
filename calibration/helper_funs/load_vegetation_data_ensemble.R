#----------------------------------------------------------------------------------------------------------------------
#The vegetation data by Paciorek et al. consists of 250 draws from the posterior distribution of a CAR model
# we want to draw samples from this distribution (they have to be internally consistent), i.e. the same
# sample for all sites

#draw 100 samples from the posterior and use thos for calibration
library(ncdf4)
veg <- nc_open('~/workflow_stepps_calibration/calibration/data/SetTreeComp_Level2_v1.0.nc')
veg.names <- names(veg$var) #gives the names of the genera available 
#load a taxon at a time
veg.comp <- lapply(veg.names,function(x) {
  ncvar_get(veg,varid=x)
})

names(veg.comp) <- veg.names

#find coordinates
meters.east <- veg$dim$x$vals
meters.north <- veg$dim$y$vals




#----------------------------------------------------------------------------------------------------------------------
#prepare that data fro STEPPS
#Simon uses structure x coordinate y coordinate,region, water, and then vegetation

#first let's prepare


sample.nr <- sample(1:250,100,replace=FALSE)
#----------------------------------------------------------------------------------------------------------------------


total.dat <- lapply(sample.nr,function(x) {
  veg.post.sample <- sapply(names(veg.comp),function(z) {
    veg.comp[[z]][,,x]
  })

total.dat <- as.data.frame(cbind(rep(meters.east,length(meters.north)),rep(meters.north,each = length(meters.east)),veg.post.sample))
colnames(total.dat)[1:2] <-c('meters.east','meters.north')
total.dat
})
# sapply(names(veg.comp),function(x){
#   image.plot(meters.east,meters.north,matrix(ncol = length(meters.north),total.dat[[x]]),zlim=c(0,1),main=x)
# })
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
write.csv(sample.nr,'~/vegetation_ensemble/sample_number_vegetation_ensemble_neus.csv')
#saveRDS(total.dat,'~/workflow_stepps_calibration/calibration/data/vegetation_samples.RDS')
