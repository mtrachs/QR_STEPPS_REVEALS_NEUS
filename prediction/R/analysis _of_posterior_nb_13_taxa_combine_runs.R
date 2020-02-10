library(rstan)
library(RColorBrewer)
library(fields)
library(sp)
library(maptools)
library(colorRamps)
#-----------------------------------------------------------------------------------------------
setwd(paste0(wd,'prediction/'))
help.fun.loc <- 'utils/'
data.loc <- 'output_nb/'
plot.loc <- 'plots_nb/'
#-----------------------------------------------------------------------------------------------
us.shp <- readShapeLines('data/map_data/us_alb.shp',proj4string=CRS('+init=epsg:3175'))
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Oak','Hemlock'))
K <- length(taxa)

prediction.files <- list.files(data.loc)

for(x in prediction.files) { 
  print(x)
  fit <- read_stan_csv(paste(data.loc,x,sep=''))
  post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  #dim(post)
  var_names<- colnames(post[,1,])
  par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
  idx_complete <- rowSums(post)>0
  post <- array(post[idx_complete,,],dim=c(sum(idx_complete),1,dim(post)[3]))
  colnames(post[,1,]) <- var_names
  gauss.process <- grep('g',par_names)
  gauss.process <- gauss.process[-1]
  
  #---------------------------------------------------------------------------------------------------------------------
  post.exp.interim <- exp(post[251:(dim(post)[1]-1),1,gauss.process])
  if(x==prediction.files[1]){ 
    post.exp <- post.exp.interim
  }else{
    post.exp <- rbind(post.exp,post.exp.interim)
  }
  
}
  ash.index <- seq(1,ncol(post.exp),K)
  
  summary.proportion <- 
    sapply(ash.index,function(z){
      exp.site <- cbind(post.exp[,z:(z+(K-1))])
      prop.site <- exp.site/rowSums(exp.site)
      prop.mean.site <- colMeans(prop.site)
      uncertainty <- apply(prop.site,2,function(x) quantile(x,probs=c(0.05,0.5,0.95)))
      data.frame(mean.proportion = prop.mean.site,median = uncertainty['50%',],lb = uncertainty['5%',],ub = uncertainty['95%',])
    })
  summary.proportion <- t(summary.proportion)
  
  mean.proportion <- matrix(ncol=K,unlist(summary.proportion[,'mean.proportion']),byrow=TRUE)
  colnames(mean.proportion) <- taxa
  
  
  
  median.proportion <- matrix(ncol=K,unlist(summary.proportion[,'median']),byrow=TRUE)
  colnames(median.proportion) <- taxa
  
  r <- median.proportion
  saveRDS(r,paste('median_prediction',x,'.RDS',sep=''))
  # 
  lb.proportion <- matrix(ncol=K,unlist(summary.proportion[,'lb']),byrow=TRUE)
  colnames(lb.proportion) <- taxa
  
  ub.proportion <- matrix(ncol=K,unlist(summary.proportion[,'ub']),byrow=TRUE)
  colnames(ub.proportion) <- taxa
  
  uncertainty <- ub.proportion - lb.proportion
  
  
  figure.name <- strsplit(x,'.csv')[[1]][1]
  
  
  load(paste0(wd,'calibration/data/for_prediction_110_sites_only_abies_new_species.RData')) 
  sputm <- SpatialPoints(pollen_coords, proj4string=CRS("+init=epsg:4326"))
  spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
  pollen_coord_us <- spgeo@coords
    
  #x1 <- strsplit(x,'.csv')[[1]][1]
  
  pdf(paste(plot.loc,'combined_110_umw.pdf',sep=''),width = 10,height = 5) #somewhat strange filename, does this refer to spatial parameters?
  par(mfrow=c(1,2),oma=c(1,1,1,2),cex = 0.8,cex.axis = 0.8)
  sapply(taxa, function(y){
    
    
    breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
    categories <- cut(mean.proportion[,y],breaks,include.lowest = TRUE,labels = FALSE)
    #cat.sort <- sort(levels(categories))
    #categories <- matrix(categories,ncol=ncol(mean.proportion))
    colours <- rev(brewer.pal(10,'RdYlBu'))
    #colours <- factor(1:10)#factor(categories)
    #colours <- tim.colors(length(breaks))
    colours.plot <- colours[categories]
    
    east <- sort(unique(veg_coords[,'meters.east']))
    north <- sort(unique(veg_coords[,'meters.north']))
    
    breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
    breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)
    
    plot(us.shp,xlim=range(east),ylim=range(north),main = paste(y,'STEPPS'))
    image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
               breaks= breaks1,lab.breaks = breaks2,cex.axis = 0.8,add=TRUE)
    points(veg_coords,col=colours.plot,pch = 15)
    points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
    plot(us.shp,add=TRUE)
    #})
    #----------------------------------------------------------------------------------------------------------------------
    #veg_agg_matrix <- as.matrix(veg_agg)
    veg <- as.matrix(r)
    
    categories <- cut(veg[,y],breaks,include.lowest = TRUE,labels = FALSE)
    
    colours.plot <- colours[categories]
    
    
    #sapply(taxa, function(x){
    plot(us.shp,xlim=range(east),ylim=range(north),main = paste(y,'Paciorek et al. (2016)'))  
    image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
               breaks= breaks1,lab.breaks = breaks2,cex.axis = 0.8,add = TRUE)
    points(veg_coords,col=colours.plot,pch = 15)
    points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
    plot(us.shp,add=TRUE)
  })
  dev.off()
#})
  
data.loc <- 'process_output/'  
saveRDS(mean.proportion,paste(data.loc,'combined_mean.RDS',sep=''))  
saveRDS(median.proportion,paste(data.loc,'combined_median.RDS',sep=''))
saveRDS(lb.proportion,paste(data.loc,'combined_lb.RDS',sep=''))
saveRDS(ub.proportion,paste(data.loc,'combined_ub.RDS',sep=''))


