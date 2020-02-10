library(stepps)#
library(dplyr)#
library(DT)#
library(neotoma)#
library(sp)#
library(fields)#
library(readr)
library(rioja)

###########################################################################################
# setwd('~/stepps-cal/R/')
# files <- list.files()
# sapply(files,function(x) source(x))
#-------------------------------------------------------------------------------------------------------------------
setwd(paste0(wd,'calibration/'))
help.fun.loc <- 'helper_funs/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-------------------------------------------------------------------------------------------------------------------

#load vegetation data using other code

if('veg_mean1.RDS'%in%list.files(data.loc)){
  veg_mean <- readRDS(paste(data.loc,'veg_mean1.RDS',sep=''))
  coords.neus <- matrix(ncol=2,c(rep(c(-79,-67),each =2),rep(c(39.5,47.5),2)))
  coords.neus <- as.data.frame(coords.neus)
  colnames(coords.neus) <- c('Lon','Lat')
  
  pol_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                       '+init=epsg:4326', 
                       '+init=epsg:4326')
  
  veg_box <- bbox_tran(coords.neus, '~ Lon + Lat',
                       '+init=epsg:4326', 
                       '+init=epsg:3175')
  
  reconst_grid <- build_grid(veg_box, resolution = 8000, proj = '+init=epsg:3175')

  
  }else {
  source(paste(help.fun.loc,'load_vegetation_data_median.R',sep=''))
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
  reconst_grid <- build_grid(veg_box, resolution = 8000, proj = '+init=epsg:3175')

  setwd(paste0(wd,'calibration/'))
  saveRDS(veg_mean,paste(data.loc,'veg_mean1.RDS',sep=''))
}





source(paste(help.fun.loc, 'get_meta_data_cal.R',sep=''))#neotoma::get_dataset(loc = pol_box, datasettype = 'pollen')#perhaps use Andrias Code for that
datasets <- meta 

#--------------------------------------------------
#load site.ids of sites that have data between present and 4000 cal BP

site_ids_us <- read.table(paste0(wd,'expert_elicitation/data/site_ids.txt'),header=TRUE)
ind.plot <- read.table(paste0(wd,'expert_elicitation/data/plot_index.txt'))

# find site ids that were plotted and ultimately used
site_ids_us <- site_ids_us[[1]][ind.plot[[1]]]
site_ids_us <- c(site_ids_us[1:70],4557,site_ids_us[71:122])

#----------------------------------------------------------------------------------------------------------

if(!'downloads.rds' %in% list.files('data/')) {
  downloads <- neotoma::get_download(datasets)
  saveRDS(downloads, paste(data.loc,'downloads.rds',sep=''))
} else {
  downloads <- readRDS(paste(data.loc,'downloads.rds',sep=''))
}


#find data that has been plotted
downloads.clean <-
  lapply(1:length(downloads), function(x) {
    daten <- NA
    if(downloads[[x]]$dataset$dataset.meta$dataset.id %in% site_ids_us) daten <- downloads[[x]]
    daten
  })

sites.pull.use <- which(sapply(1:length(downloads.clean), function(x) unique(!is.na(downloads.clean[[x]]))))

downloads.clean <- lapply(sites.pull.use, function(x) downloads.clean[[x]])

#try to sort site ids so that sample.id is ascending (is alrady done in fact not necessary)
site.ids.downloads.clean <- sapply(1:length(downloads.clean), function(x) 
  downloads.clean[[x]]$dataset$dataset.meta$dataset.id)


#site.ids.downloads.clean[order(site.ids.downloads.clean)]
downloads.clean1 <- lapply((1:length(downloads.clean))[order(site.ids.downloads.clean)],function(x) downloads.clean[[x]])
#-----------------------------------------------------------------------------------------------------------------------
#load evaluation of elecitation exercise
source(paste(help.fun.loc,'evaluate_elicitation_certainty.R',sep=''))

#(1:length(sites.use))[sites.use]
#sample.use


#--------------------------------------------------------------------------------------------------------------------
# find sites that are outside the domain of vegetation
#immediately remove the three sites that were not used in elicitation
if('meta_data.RDS'%in%list.files(paste0(wd,'expert_elicitation/data/'))){
  meta.data.neus <- readRDS(paste0(wd,'expert_elicitation/data/meta_data.RDS'))
}else {
source(paste0(wd,'expert_elicitation/elicitation_helper_funs/get_meta_data.R'))
  saveRDS(meta.data.neus,paste0(wd,'expert_elicitation/data/meta_data.RDS'))
}  
  
ind.plot <- ind.plot[[1]]
site.real <- meta.data.neus$site[order(meta.data.neus$datasetID)][ind.plot]
long.real <- meta.data.neus$long[order(meta.data.neus$datasetID)][ind.plot]
lat.real <- meta.data.neus$lat[order(meta.data.neus$datasetID)][ind.plot]


#transform to us coordinates
sputm <- SpatialPoints(cbind(long.real,lat.real), proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))

lakes.coord <- spgeo@coords
#calcualte distances
dist.lakes <- paldist2(lakes.coord,cbind(veg_mean$meters.east,veg_mean$meters.north), dist.method="euclidean")
min.dist.lakes <- apply(dist.lakes,1,min)
min.dist.lakes.index <- which(min.dist.lakes<8000)
#----------------------------------------------------------------------------------------------------------------------
#total removal criterion
total.exclude <- ((min.dist.lakes > 8000))#|((long.real> (-70.7))& (lat.real<43))|((long.real> (-70))& (lat.real<44.5))|
  # ((long.real> (-68)))|((long.real < (-69.5))& (lat.real>45.25))|((lat.real < 39.8))|
   #((long.real < (-75))& (lat.real>44.25)))

map('state',xlim=c(-82,-67),ylim=c(38,50))
points(cbind(long.real,lat.real)[total.exclude==FALSE,],pch = 16)
points(cbind(long.real,lat.real)[total.exclude,],pch = 16, col = 3)

#which sites to use
coordinates.use.index <- which(total.exclude==FALSE)
#------------------------------------------------------------------------------------------------------------------
#first remove sites that were not plotted from data.NE

#take intersect between sites chosen by experts and sites that have good coordinates
index.coordinate.expert  <- intersect(coordinates.use.index,site.sample[,'site.index'])
#object site.sample comes from evaluate_elicitation_certainty.R file  

#break down site index and sample number so that only places with good coordinates are used
site.sample.good.coordinate <- as.data.frame(site.sample[site.sample[,'site.index']%in%index.coordinate.expert,]) 

test <- 
  sapply(site.sample.good.coordinate$site.index,function(x) {#site.sample.good.coordinate$site.index
    sam.index <- site.sample.good.coordinate$sample_number[site.sample.good.coordinate$site.index==x] 
    counts <- downloads.clean1[[x]]$counts[sam.index,]
    depth <-  downloads.clean1[[x]]$sample.meta$depth[sam.index]
    site.id <- downloads.clean1[[x]]$dataset$site$site.id
    dataset.id <- downloads.clean1[[x]]$dataset$dataset.meta$dataset.id
    site.name <- downloads.clean1[[x]]$dataset$site$site.name
    lon <- downloads.clean1[[x]]$dataset$site$long
    lat <- downloads.clean1[[x]]$dataset$site$lat
    complete.data <- data.frame(site.id = as.character(site.id),
                                dataset.id = dataset.id,
                                site.name = as.character(site.name),
                                lon = lon,
                                lat = lat,
                                depth = depth,
                                sample.number = sam.index,
                                t(counts))
  })


#merge taxa for each site??


########################################################################################################################
#load samples from two additional lakes ## in fact four additional sites 
########################################################################################################################
additional_dataset_ids <- c(15352,15598,2289,523)
add_samples <- c(8,9,2,3)

two.sites <- lapply(additional_dataset_ids,function(x) {# run over the site indexes found above
  # download data for a specified site
  data.download <- get_download(x) 
  # get depth of a sample
  depth <- data.download[[1]]$sample.meta$depth
  lon <- data.download[[1]]$dataset$site$long
  lat <- data.download[[1]]$dataset$site$lat 
  name.site <- data.download[[1]]$dataset$site$site.name
  site.id <- data.download[[1]]$dataset$site.data$site.id
  #find chronologies of this core
  #some cores have more than one chronology stored, we have to check if any of these chronologies contains a date in the 
  #period we are interested in 
  
  #load all counts in the first place
  all.counts <- data.download[[1]]$counts
  #find Lycopdoium, Eucalyptus and Microspheres
  ind.lyco <- grep('Lycopodium',colnames(all.counts))
  ind.euca <- grep('Eucalyptus',colnames(all.counts))
  ind.micro <- grep('Micros',colnames(all.counts))
  ind.tot <- c(ind.lyco,ind.euca,ind.micro)
  if(length(ind.tot)>0) all.counts <- all.counts[,c(-ind.tot)]
  counts <- as.matrix(all.counts)
  # if we only have on sample this is written in a n x 1 matrix, I want to change that
  if(min(dim(counts))==1) counts <- t(counts)
  # turn counts in to percentages
  if(nrow(counts)>1) percentages <- round(100 * counts/rowSums(counts),2)
  # turn counts in to percentages
  if(nrow(counts)==1) percentages <- round(100 * counts/sum(counts),2)
  # exclusion rule we need at least 2% in one entry
  include.clean <- apply(percentages,2,function(x) max(x)>2) 
  # applyexclusion rule
  percentages.clean <- percentages[,include.clean]
  #extract depths used later on (within a certain chronology)  
  # load chornology of selected data 
  #else it also takes depths with NA
  #return data for later use
  list(counts = counts, percentages = percentages,percentages.clean = percentages.clean,
       dataset_idx = x,site_idx = site.id,depth = depth,lat = lat, lon=lon,site.name = name.site)
})

l.data <- length(test)

for(i in 1:length(add_samples)) {
  sample.number <- add_samples[i]
  daten <- two.sites[[i]]
  test[[l.data+i]] <- data.frame(site.id = as.character(daten$site_idx),
                                 dataset.id = daten$dataset_idx,
                                 site.name = as.character(daten$site.name),
                                 lon = daten$lon,
                                 lat = daten$lat,
                                 depth = daten$depth[sample.number],
                                 sample.number = sample.number,
                                 t(daten$counts[sample.number,]))
 
}

########################################################################################################################
# add data not in neotoma
########################################################################################################################

l.data <- length(test)

data.not.in.neotoma <- readRDS(paste(data.loc,'settlement_era_data_not_in_neotoma.RDS',sep=''))
coord.New.England <- readRDS(paste(data.loc,'coordinates_New_England.RDS',sep=''))
coord.not.in.neotoma.NEUS <- readRDS(paste(data.loc,'coordinates_not_in_neotoma_NEUS.RDS',sep=''))
colnames(coord.not.in.neotoma.NEUS) <- c('lon','lat')

coord.additional.data <- rbind(coord.not.in.neotoma.NEUS,coord.New.England)



for(i in 1:nrow(coord.additional.data)) {
  #sample.number <- add_samples[i]
  #daten <- two.sites[[i]]
  test[[l.data+i]] <- data.frame(site.id = NA,
                                 dataset.id = NA,
                                 site.name = names(data.not.in.neotoma)[i],
                                 lon = coord.additional.data [i,'lon'],
                                 lat = coord.additional.data [i,'lat'],
                                 depth = NA,
                                 sample.number = NA,
                                 data.not.in.neotoma[[i]])
  
}



#pollen.test <- as.data.frame(t(test[[1]]))
pollen.test <- as.data.frame(test[[1]])


for (i in 2:length(test)) {
  pollen.test <- analogue::join(as.data.frame(pollen.test),as.data.frame(test[[i]]),split=FALSE) # gives a few warnings they are ok
}



##################################################################################################################
#remove sample for Goose Bay Marsh, this sample does only contain 67 pollen grains while 166 grains of pinaceae undiff are removed
gbm <- which(pollen.test$dataset.id==15732)
#pollen.test[gbm,]
pollen.test <- pollen.test[-gbm,]
#################################################################################################################
#write.csv(pollen.test[,c('site.name','site.id','dataset.id','lon','lat','depth','sample.number')],
#          paste(data.loc,'meta_data_calibration_samples.csv',sep=''),row.names = FALSE)




write.csv(sort(colnames(pollen.test)),'~/workflow_stepps_calibration/calibration/data/taxon_names.csv')
  

#make taxa translation table 
source(paste(help.fun.loc,'taxa_translation_only_abies.R',sep=''))
pol_table <- readr::read_csv(paste(data.loc,'taxon_names_translated.csv',sep=''))


#perhaps need to comment this
calib_trans <- translate_taxa(pollen.test, 
                              pol_table,
                              id_cols = colnames(pollen.test)[1:7])


#save data to be able to add sample to age information for age-depth modelling
saveRDS(calib_trans,paste(data.loc,'settlement_era_pollen.RDS',sep=''))


#-------------------------------------------------------------------------------------------------------------------
veg_table <- readr::read_csv(paste(data.loc,'/veg_trans_edited_only_abies.csv',sep=''))
veg_trans <- translate_taxa(veg_mean, veg_table ,id_cols = colnames(veg_mean)[1:2])


#-------------------------------------------------------------------------------------------------------------------
veg_table <- to_stepps_shape(veg_trans,   '~ meters.east + meters.north',      '+init=epsg:3175')
pol_table <- to_stepps_shape(calib_trans, '~ lon + lat', '+init=epsg:4326')

target_taxa <- colnames(calib_trans)[8:ncol(calib_trans)]

source(paste(help.fun.loc,'prep_input_modified.R',sep=''))
test_neus <- prep_input(veg = veg_table, 
                 pollen = pol_table, 
                 target_taxa = target_taxa,
                 grid   = reconst_grid,hood = 7e+05)

test_neus$d <-  round(test_neus$d/1e+06,5)
test_neus$d_pot[,1] <- test_neus$d_pot[,1]/1e+06 

N_cores <- test_neus$N_cores
N_cells <- test_neus$N_cells
y <- test_neus$y
r <- test_neus$r
d <- test_neus$d
idx_cores <- test_neus$idx_cores
idx_hood <- test_neus$idx_hood
d_pot <- test_neus$d_pot
N_pot <- test_neus$N_pot
N_hood <- test_neus$N_hood
num_sites <- dim(pol_table)[1]



oh_new_pollen <- rowSums(y[,c("Other hardwood","Elm","Poplar")])
y[,"Other hardwood"] <- oh_new_pollen
y <- y[,(colnames(y)%in%c("Elm","Poplar"))==FALSE]
K <- ncol(y)

oh_veg <- rowSums(r[,c("Other hardwood","Elm","Poplar")])
r[,"Other hardwood"] <- oh_veg
r <- r[,(colnames(r)%in%c("Elm","Poplar"))==FALSE]
r <- r/rowSums(r)

library(rstan)
stan_rdump(list = names(test_neus), file = paste(data.loc,'/for_prediction_',num_sites,'_sites_only_abies_new_species.dump',sep=''))

veg_coords <- veg_table@coords
pollen_coords <- pol_table@coords
save(list = c(names(test_neus),'veg_coords','pollen_coords'), 
     file = paste(data.loc,'for_prediction_',num_sites,'_sites_only_abies_new_species.RData',sep=''))
#-------------------------------------------------------------------------------------------------------------------
