library(ncdf4)
veg <- nc_open(paste0(data.loc,'SetTreeComp_Level2_v1.0.nc'))
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
    apply(veg.comp[[z]],c(1,2),function(x) median(x,na.rm=TRUE))
  })

total.dat <- as.data.frame(cbind(rep(meters.east,length(meters.north)),rep(meters.north,each = length(meters.east)),
                                 veg.mean.post))
colnames(total.dat)[1:2] <-c('meters.east','meters.north')

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# 
#sp::coordinates(total.dat) <- ~ meters.east + meters.north

# pol_box <- bbox_tran(total.dat, '~ meters.east + meters.north',
#                      '+init=epsg:3175', 
#                      '+init=epsg:4326')
# 
# veg_box <- bbox_tran(total.dat, '~ meters.east + meters.north',
#                      '+init=epsg:3175', 
#                      '+init=epsg:3175')
# 
# datasets <- neotoma::get_dataset(loc = pol_box, datasettype = 'pollen',ageold = )
# 
# 
# 
# #-----------------------------------------------------------------------------------------------------------------------
# #try to find projection of 
# library(rgdal)
# library(maps)
# coords.neus <- matrix(ncol=2,c(rep(c(-80.5,-67),each =2),rep(c(39.5,47.5),2)))
# 
# spgeo <- SpatialPoints(coords.neus, proj4string=CRS('+init=epsg:4326'))
#                        #("+proj=utm +zone=15 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# spalbers <- spTransform(spgeo, CRS('+init=epsg:3174'))
# 
# coord.albers <- coordinates(spalbers)





# library(fields)
# 
# breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# colours <-  tim.colors(length(breaks)-1)
# 
# pdf("C://Work/Paleon/Datasets/vegetation_NEUS1.pdf",width = 11.69,height = 8.29)
#   sapply(veg.names,function(x) {
#     taxon <- matrix(ncol=length(meters.north),nrow=length(meters.east)) 
#   
# 
#       for (i in 1: length(meters.east)) {
#         for(ii in 1:length(meters.north)) {
#           taxon[i,ii] <- mean(veg.comp[[x]][i,ii,],na.rm=TRUE) 
#         }
#      } 
# 
#     image.plot(meters.east,meters.north,taxon,main =x,xlim=c(1200000,2300000),ylim=c(350000,1300000),
#                breaks = breaks,col=colours)
#   })
# dev.off()
# 
# 
# 
# #check how data goes into a matrix
# pdf("C://Work/Paleon/Datasets/vegetation_NEUS2.pdf",width = 11.69,height = 8.29)
# sapply(veg.names,function(x) {
#   taxon <- matrix(ncol=length(meters.north),nrow=length(meters.east),veg.mean.post[,x]) 
#     image.plot(meters.east,meters.north,taxon,main =x)
# })
# dev.off()
# 
# 
# 
# 
# #----------------------------------------------------------------------------------------------------------------------
# #assign colours
# 
# breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# 
# categories <- cut(veg_table$Oak,breaks,include.lowest=TRUE)
# colours <-  tim.colors(length(breaks)-1)
# 
# colour <- vector(length=length(categories))
# 
# for (i in 1:length(unique(categories))) {
#   colour[categories==sort(unique(categories))[i]] <- colours[i]
# }
# 
# colour <- replace(colour,colour==FALSE,'white') 
#   
# plot(veg_table@coords[,'meters.east'],veg_table@coords[,'meters.north'],col=colour,pch = 15)
# 
# 
# 
# 
# #--------------------------------------------------------------------------------------------------------------------
# breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# 
# categories <- cut(veg_mean$Oak,breaks,include.lowest=TRUE)
# colours <-  tim.colors(length(breaks)-1)
# 
# colour <- vector(length=length(categories))
# 
# for (i in 1:length(unique(categories))) {
#   colour[categories==sort(unique(categories))[i]] <- colours[i]
# }
# 
# colour <- replace(colour,colour==FALSE,'white') 
# 
# plot(veg_mean$meters.east,veg_mean$meters.north,col=colour,pch = 15)
# 
# 
# 
# 
# #---------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------
# breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# 
# categories <- cut(total.dat$Oak,breaks,include.lowest=TRUE)
# colours <-  tim.colors(length(breaks)-1)
# 
# colour <- vector(length=length(categories))
# 
# for (i in 1:length(unique(categories))) {
#   colour[categories==unique(categories)[i]] <- colours[i]
# }
# 
# colour <- replace(colour,colour==FALSE,'white') 
# 
# plot(total.dat$meters.east,total.dat$meters.north,col=colour,pch = 15)
# 
# 
# #----------------------------------------------------------------------------------------------------------------------
# #check what happens with lapply
# veg.mean.post.lapply <- 
#   lapply(names(veg.comp),function(z) {
#     apply(veg.comp[[z]],c(1,2),function(x) round(mean(x),3))
#   })
# 
# 
# names(veg.mean.post.lapply) <- veg.names
# 
# image.plot(meters.east,meters.north,veg.mean.post.lapply$Oak)#looks good
# 
# #put data into matrix manually
# oak.test <- matrix(veg.mean.post.lapply$Oak,ncol=1)
# 
# oak.test[1:1000]
# 
# min(which(!is.na(oak.test)))
# 
# veg.mean.post.lapply$Oak[,2]#meters east change while meters north remain constant
# 
# #-------
# breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
# 
# categories <- cut(oak.test,breaks,include.lowest=TRUE)
# colours <-  tim.colors(length(breaks)-1)
# 
# colour <- vector(length=length(categories))
# 
# for (i in 1:length(unique(categories))) {
#   colour[categories==unique(categories)[i]] <- colours[i]
# }
# 
# colour <- replace(colour,colour==FALSE,'white') 
# 
# plot(total.dat$meters.east,total.dat$meters.north,col=colour,pch = 15)
# 
