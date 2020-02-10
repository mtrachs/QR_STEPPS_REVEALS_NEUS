library(raster)
library(sp)
library(DISQOVER)
library(reshape2)
library(ggplot2)
library(dplyr)
library(neotoma)
library(maps)
library(RColorBrewer)
library(fields)

#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)

# # load package
# load('DISQOVER/R/sysdata.rda')

setwd(paste0(wd,'/Reveals_NEUS/'))
data.loc <- paste0(wd,'calibration/data/') 
help.fun.loc <- paste0(wd,'calibration/helper_funs/')

ena <- TRUE

##################################################################################################################################################
## pull pollen data for NA
##################################################################################################################################################

# sum samples within a time bin for each site 

# make grid for NA (or ENA)
source('R/make_grid.R')
#-------------------------------------------------------------------------------------------------------------------
#find coordinates of NEUS in Degrees
#-------------------------------------------------------------------------------------------------------------------
# load vegetation data prepared for STEPPS
load(paste(data.loc,'elicitation_neus_certainty_median_110_sites_only_abies_new_species.RData',sep=''))
pollen_coords <- as.data.frame(pollen_coords)
pollen_bin <- data.frame(pollen=y,lon = pollen_coords$lon,lat = pollen_coords$lat)



#-----------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#lake sizes
#-----------------------------------------------------------------------------------------------------------------
meta_dat <- as.data.frame(read.csv(paste0(data.loc,'pollen_elicitation_NEUS.csv')))

site.names <- meta_dat$site.name#readRDS(paste(data.loc,'lake_names.RDS',sep=''))
site.names <- as.character(site.names)
site.ids <- meta_dat$site.id #readRDS(paste(data.loc,'site_ids.RDS',sep=''))

site.names[site.names=='black'] = 'Black not neotoma'
site.names[site.names=='North Round Pond'] = 'North Round Pon.'

lake_size <- read.csv(paste0(wd,'Reveals_NEUS/data/assigned_areas_1.6.csv'))
# would be better to take a 
lake_index_NEUS <- lake_size$stid%in%site.ids
lake_sizes_NEUS <- lake_size[lake_index_NEUS,]



######################################################################
#find sites in neotoma without sizes
######################################################################
lake_index_NEUS_size_na <- !(site.ids%in%lake_size$stid)
site_names_NEUS_na <- site.names[intersect(which(lake_index_NEUS_size_na),which(!is.na(site.ids)))]
#######################################################################
map('state')
points(lake_sizes_NEUS$long,lake_sizes_NEUS$lat,pch = 16,cex = 0.5)


#lake_sizes_NEUS <- lake_sizes_NEUS[((lake_sizes_NEUS$long>(-90))&(lake_sizes_NEUS$lat>38)),]
#lake_size <- lake_sizes_NEUS$area

###################################################################################################################
#read lake sizes given in Oswald et al. 
###################################################################################################################
lake_sizes_new_england <- read.csv(paste0(wd,'new_pollen_data/CSV_new/lake_sizes_Oswald_et_al.csv'))
lake_size_estimated <- read.csv(paste0(wd,'new_pollen_data/CSV_new/additional_estimated_sizes.csv'))
lake_sizes_new_england1 <- lake_sizes_new_england[,c('Lake.Name','size')]
lake_sizes_new_england2 <- rbind(lake_sizes_new_england1,lake_size_estimated)
lake_sizes_new_england2$Lake.Name <- as.character(lake_sizes_new_england2$Lake.Name)
lake_sizes_new_england2$Lake.Name[as.character(lake_sizes_new_england2$Lake.Name)=="Black"] <- 'Black not neotoma'
lake_sizes_new_england2$Lake.Name[as.character(lake_sizes_new_england2$Lake.Name)=="North Round Pond"] <- 'North Round Pon.'


# some of the sites by Oswald et al. are already in NEeotoma, others aren't. 

#find sites not in neotoma that have sizes
site.names_not_neotoma <- site.names[intersect(which(lake_index_NEUS_size_na),which(is.na(site.ids)))]

#find additional sites with lake sizes
index.not.neotoma <- unlist(sapply(site.names_not_neotoma,function(x) { 
  grep(x,lake_sizes_new_england2$Lake.Name,ignore.case=TRUE)
  })) # this index is an index within the new England lake size dataset not the sites...

lake_sizes_new_england_not_neotoma <- lake_sizes_new_england2[sort(index.not.neotoma),]
lake_sizes_new_england_surplus <- lake_sizes_new_england2[-c(sort(index.not.neotoma)),]

#################################################################################################################################
#was used to find lakes without lake sizes 
#find the index within the lake names
# index.lake.names <- unlist(sapply(lake_sizes_new_england_not_neotoma$Lake.Name,function(x) { 
#   grep(x,site.names_not_neotoma,ignore.case=TRUE)
# })) # this index is an index within the new England lake size dataset not the sites...
# 
# 
# sites.not.neotoma.missing <- site.names_not_neotoma[!(1:length(site.names_not_neotoma))%in%index.lake.names]
# 
# #look for missing sites in neotoma 
# hae  <- unlist(sapply(sites.not.neotoma.missing,function(x) { 
#   grep(x,lake_size$name,ignore.case=TRUE)
# }))

#lake_size$name[hae]
###################################################################################################################################
vector_not_neotoma <- unlist(sapply(1:length(lake_sizes_new_england_not_neotoma$Lake.Name),function(x) { 
  grep(lake_sizes_new_england_not_neotoma$Lake.Name[x],site.names,ignore.case=TRUE)
}))

vector_neotoma <- unlist(sapply(1:length(lake_sizes_NEUS$stid),function(x) { 
  ll <-which(lake_sizes_NEUS$stid[x]==site.ids)
  if((length(ll)>1)&(x%%2==0)){
    ll <- ll[1]
  }
  if((length(ll)>1)&(x%%2==1)){
    ll <- ll[2]
  }
  ll
}))



lake.sizes.final <- as.data.frame(matrix(ncol=2,nrow=length(site.names)))
colnames(lake.sizes.final) <- c('Lake.Name','area')

for(i in 1:length(vector_not_neotoma)){
  lake.sizes.final[vector_not_neotoma[i],1] <- as.character(lake_sizes_new_england_not_neotoma$Lake.Name[i])
  lake.sizes.final[vector_not_neotoma[i],2] <- lake_sizes_new_england_not_neotoma$size[i]
}


for(i in 1:length(vector_neotoma)){
  lake.sizes.final[vector_neotoma[i],1] <- as.character(lake_sizes_NEUS$name[i])
  lake.sizes.final[vector_neotoma[i],2] <- lake_sizes_NEUS$area[i]
}

cbind(site.names,lake.sizes.final)
write.csv(cbind(site.ids,lake.sizes.final),paste0(data.loc,'lake_sizes_NEUS.csv'))

#global.index <- #(83:110)[sort(index.lake.names)] #this is probably wrong
#neotoma.index <- site.names%in%lake_sizes_NEUS$site.name

# see if we can find sizes for sites that are not in neotoma for instance from Jackson and Whitehead

# 




# check if there is additional data 
#site_names_NEUS_na%in%lake_sizes_new_england_surplus$Lake.Name # no additional lake sizes
##-----------------------------------------------------------------------------------------------------------------
#function to generate random lake size
random_lake_size <- function(lake_size){
  library(truncnorm)
  lake_radius = rtruncnorm(1,
                           a=0.5,
                           mean=median(lake_size, na.rm=TRUE),
                           sd=IQR(lake_size, na.rm=TRUE))
  
  return(lake_radius)
}

#---------------------------------------------------------------------------------------------------------------
#produce random lake sizes
lakes_no_size <- sum(is.na(lake.sizes.final$area))
rand_la_size <- replicate(lakes_no_size,random_lake_size(lake_size = lake.sizes.final$area))

lake.sizes.final$area[is.na(lake.sizes.final$area)] <- rand_la_size 

#sort_lake_names <-  site.names[order(site.names)]
#sort_lake_sizes <- lake_sizes_NEUS[order(lake_sizes_NEUS$site.name),]

#a few sites have two pollen samples
#total_lake_sizes <- vector(length = length(site.names))
#total_lake_sizes[global.index] <- lake_sizes_new_england_not_neotoma$size 

#for (i in 1:length(sort_lake_sizes$site.name)) {
#  total_lake_sizes[as.character(sort_lake_names)==as.character(sort_lake_sizes$site.name[i])]   <- sort_lake_sizes$area[i] 
#}

#total_lake_sizes[total_lake_sizes==0] <- replicate(sum(total_lake_sizes==0),random_lake_size(lake_size = lake_size))


# in HA; convert to radius in m
# pi * r * r
lake_r = sqrt(lake.sizes.final$area*0.01 / pi)*1000
lake_radius <- lake_r
lake_radius[lake_radius<10] <- 10
#------------------------------------------------------------------------------------------------------------------
#have to find the data that we want to use at the end of the day!!!
#data used for calibration plus data used for 
#function that evaluates sites
help.fun.loc <- paste0(wd,'calibration/helper_funs/')
#source(paste(help.fun.loc,'evaluate_elicitation_certainty.R',sep=''))
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------







for(reg.cutoff in c(50000,400000,100000)) { 
  for(larix in c('larch','no_larch')) {
    for(la_size in c('correct','random')) {

taxa <- colnames(y)
taxa[taxa%in%c('Tamarack','Other conifer','Other hardwood')] <-c('Fir','Alder','Larch')
taxa <- taxa[!(taxa%in%'Chestnut')]
if(larix=='no_larch') taxa <- taxa[!(taxa%in%'Larch')]

grid <- make_grid(pollen_bin, coord_fun = ~ lon + lat, projection = '+init=epsg:4326', resolution = 1)

cell_id <- raster::extract(grid, pollen_bin[,c('lon', 'lat')])

grid_coords <- rasterToPoints(grid)
# 
if(la_size=='random') la_radius <- sample(lake_radius,length(lake_radius),replace=FALSE)
if(la_size!='random') la_radius <- lake_radius

#look at this...
pollen_bin1 <- data.frame(la_radius, pollen_bin)


##################################################################################################################################################
## read in and prep ppes and svs
##################################################################################################################################################
ppes = readRDS(paste0(wd,'Reveals_NEUS/data/PPEs_agg.RDS'))
ppes = ppes[which(ppes$taxon %in% taxa),] #no PPE for Chestnut

#ppes[which(ppes$ppe.error == 0),'error'] = 0.01

###############################################################################################################################
## read in SVs
svs = read.csv(paste0(wd,'Reveals_NEUS/data/svs_LC6K.csv'), header=TRUE, stringsAsFactors=FALSE)
#svs$taxon[svs$taxon=='Fir'] <- 'Abies'
svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)

svs_agg = svs_agg[svs_agg$taxon %in% taxa, ]

## construct param input data frame




# sapply(c('literature','stepps'),function(ppe_estimated) {
#   sapply(c('measured','constant'),function(fall_speed) {
#    sapply(c('gpm neutral'),function(dist.w){ #,'lsm unstable'
for(ppe_estimated in 'literature') { 
  for(fall_speed in c('measured','constant')) {
    for(dist.w in c('gpm neutral')) {

    if(ppe_estimated == 'literature') {
  
      if(fall_speed=='measured') fallspeed=svs_agg$sv[svs_agg$taxon%in%ppes$taxon]
      if(fall_speed=='constant') fallspeed=rep(mean(svs_agg$sv[svs_agg$taxon%in%ppes$taxon]),length(ppes$taxon))
      reveals_inputs = data.frame(species=ppes$taxon,
                                fallspeed=fallspeed,
                                PPEs=ppes$ppe,
                                PPE.errors=ppes$error)
      reveals_inputs$species <- as.character(reveals_inputs$species)
      if(larix =='larch') reveals_inputs$species[reveals_inputs$species%in%c('Fir','Alder','Larch')] <- c('Other_hardwood','Other_conifer','Tamarack')
      if(larix =='no_larch') reveals_inputs$species[reveals_inputs$species%in%c('Fir','Alder')] <- c('Other_hardwood','Other_conifer')
      reveals_inputs <- reveals_inputs[order(reveals_inputs$species),]
      rownames(reveals_inputs) = NULL
      reveals_inputs <- reveals_inputs[order(reveals_inputs$species),]
    }

    if(ppe_estimated == 'stepps') {
      ppe_stepps <- as.data.frame(readRDS(paste0(wd,'Reveals_NEUS/data/ppe_stepps.RDS')))
      #test this line
      #ppe_stepps <-  ppe_stepps[rownames(ppe_stepps)%in%taxa,]
      #had to add elm and poplar here if we change the stepps file, this should get changed again 
      ppe_stepps <- ppe_stepps[!rownames(ppe_stepps)%in%c('Chestnut','Elm','Poplar'),]
      if(larix=='no_larch') ppe_stepps <- ppe_stepps[rownames(ppe_stepps)!='Tamarack',]
      if(larix=='larch') svs_agg$taxon[svs_agg$taxon%in%c('Alder','Fir','Larch')] <- c('Other hardwood','Other conifer','Tamarack')
      if(larix=='no_larch') svs_agg$taxon[svs_agg$taxon%in%c('Alder','Fir')] <- c('Other hardwood','Other conifer')
      #something seems odd here!!
      if(fall_speed=='measured') fallspeed=svs_agg$sv[order(svs_agg$taxon)]
      if(fall_speed=='constant') fallspeed=rep(mean(svs_agg$sv),length(taxa))#fallspeed=rep(mean(svs_agg$sv[order(svs_agg$taxon)]),length(taxa))
      reveals_inputs = data.frame(species=taxa,
                              fallspeed=fallspeed,
                              PPEs=ppe_stepps$ppe_STEPPS_mean_stand,
                              PPE.errors=ppe_stepps$ppe_STEPPS_sd_stand)
    }
    write.csv(reveals_inputs, "data/reveals_input_params.csv", row.names=FALSE)

##################################################################################################################################################
## run reveals
##################################################################################################################################################
    ids = unique(pollen_bin1$dataset)
    pol_dat = pollen_bin1


    if(larix=='larch') veg_pred <- matrix(ncol = 5*(ncol(pol_dat)-4)+2,nrow = 1) 
    if(larix=='no_larch') veg_pred <- matrix(ncol = 5*(ncol(pol_dat)-5)+2,nrow = 1) 
    colnames(veg_pred) <- NULL


    for (i in 1:nrow(pol_dat)){
  
      print(i)
      id = ids[i]
  
      counts_site = pol_dat
      colnames(counts_site)[colnames(counts_site)%in%c("pollen.Other.conifer" ,"pollen.Other.hardwood")] <- 
        c("pollen.Other_conifer" , "pollen.Other_hardwood")
      counts_site <- cbind(matrix(ncol=1, nrow(counts_site)),counts_site)
      colnames(counts_site)[1] <- 'ages'
   
      basin_radius = lake_radius
  
  
      coords_site = pollen_coords
      rownames(coords_site) = NULL
  
  
  #why do we have this issue with lake_radius.1
      if(larix=='larch') write.csv((counts_site[i,!(colnames(counts_site)%in%c('la_radius','la_radius.1','pollen.Chestnut','lon','lat'))]), 'data/reveals_input.csv', row.names=FALSE)
      if(larix=='no_larch') write.csv((counts_site[i,!(colnames(counts_site)%in%c('la_radius','la_radius.1','pollen.Chestnut','pollen.Tamarack','lon','lat'))]), 'data/reveals_input.csv', row.names=FALSE)
  
  # cycle through and estimate background veg
  # with english csv files
    a <- REVEALSinR(pollenFile = "data/reveals_input.csv",
                  pf         = "data/reveals_input_params.csv",
                  filetype   = "csv",
                  dwm        = dist.w,
                  tBasin     = "lake",
                  dBasin     = 2*round(basin_radius[i]), # diameter!
                  regionCutoff = reg.cutoff,
                  repeats      = 1000)
  
    veg = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
  
    veg$taxon = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
    veg$type  = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))
    column.names   = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), "pollen.", fixed=TRUE))[2]))
    
    veg_inter <- data.frame(coords_site[i,], t(veg$value))
    colnames(veg_inter)[3:ncol(veg_inter)] <- column.names 
  
    veg_pred = rbind(veg_pred,as.matrix(veg_inter))
    }

    #colnames(veg_pred) <- c('lon','lat',taxa[taxa!='Chestnut'])
    veg_pred <- veg_pred[-1,] 

    saveRDS(veg_pred, paste('output/unc_stand_veg_pred_ppe_',ppe_estimated,'_fallspeed_',fall_speed,
                            '_composition_',larix,'_lake_size_',la_size,'_dwm_',dist.w,'_',reg.cutoff,'_regional_vegetation.RDS',sep=''))
  
    veg_pred <- cbind(veg_pred,cell_id)
    veg_pred_agg <- aggregate(veg_pred[,-c(1:2,ncol(veg_pred))],by = list(veg_pred[,ncol(veg_pred)]),FUN = mean)
    c.names <- colnames(veg_pred_agg)
    
    coord_veg_pred <- grid_coords[grid_coords[,'layer']%in%veg_pred_agg$Group.1,]

    veg_pred_agg <- matrix(ncol = ncol(veg_pred_agg),unlist(veg_pred_agg))
    colnames(veg_pred_agg) <- c.names
    veg_pred_agg <- veg_pred_agg[,-1]
    saveRDS(veg_pred_agg, paste('output/unc_stand_veg_pred_ppe_',ppe_estimated,'_fallspeed_',fall_speed,
                            '_composition_',larix,'_lake_size_',la_size,'_dwm_',dist.w,'_',reg.cutoff,'_regional_vegetation_aggregated.RDS',sep=''))
    
    saveRDS(coord_veg_pred, paste('output/unc_stand_ratio_pred_ppe_',ppe_estimated,'_fallspeed_',fall_speed,
                                '_composition_',larix,'_lake_size_',la_size,'_dwm_',dist.w,'_',reg.cutoff,'_regional_vegetation_aggregated_coord.RDS',sep=''))
    
    #----------------------------------------------------------------------------------------------------------------------
    #look at figures
    breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
     colours <- rev(brewer.pal(10,'RdYlBu'))
     
    breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
    breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)

    east <- sort(unique(grid_coords[,'x']))
    north <- sort(unique(grid_coords[,'y']))
    #---------------------------------------------------------------------------------------------
    pdf(paste('figures/unc_stand_ppe_',ppe_estimated,'_fallspeed_',fall_speed,
              '_composition_',larix,'_lake_size_',la_size,'_dwm_',dist.w,'_reg_cutoff_',reg.cutoff,'.pdf',sep=''),height = 5,width =6)
    for (i in 1:ncol(veg_pred_agg)){
      categories <- cut(veg_pred_agg[,i]/100,breaks,include.lowest = TRUE,labels = FALSE)
      colours.plot <- colours[categories]
      par(oma=c(1,1,1,2))
      map('state',xlim=c(-81,-66.5),ylim=c(39.5,49.5))
      mtext(side = 3, lin=2.2,font = 2,text = colnames(veg_pred_agg)[i])
      image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,main = colnames(colours.plot)[i],cex.axis = 0.8,add=TRUE)
      points(coord_veg_pred,col=colours.plot,pch = 15,cex = 3)
      map('state',add=TRUE,xlim=c(-81,-66.5),ylim=c(39.5,49.5))
    }
    dev.off()
  }
}
}
}
}
}
