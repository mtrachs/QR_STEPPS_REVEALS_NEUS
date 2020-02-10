#---------------------------------------------------------------------------------------------------------------
library(neotoma) #1.7
library(analogue)#0.17-0
library(rioja)#0.9-15
library(maps)#3.2.0

#---------------------------------------------------------------------------------------------------------------
#set working directory and location of helper files
setwd(paste0(wd,'/expert_elicitation/'))
help.fun.loc <- 'helper_funs/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#---------------------------------------------------------------------------------------------------------------
# make map of the area
map(database = 'state')
abline(h = c(39.75,47.5),col=4)
abline(v = c(-85,-67),col=4)
#---------------------------------------------------------------------------------------------------------------
#this is a serious attempt to load data from the North eastern US for the time period 350 to 650 cal BP 
#---------------------------------------------------------------------------------------------------------------

#first the brackeeting ages are defined, we only consider cores that have data in this time period
older.age <- 4000
younger.age <- -70


#what is this
get_table('GeoPoliticalUnits') 

#load fossil data for NEUS (given by coordinates in loc )
new_england_fossil <- get_dataset(datasettype = 'pollen',
                                  loc=c(-85,39.75,-67,47.55)) #
#returns a long list now we are only looking for the dataset ID
#use Andria Dawson's Code to find indexes that are in New England (have to change place and name of this code) 

if('meta_data.RDS'%in%list.files(data.loc)){
  meta.data.neus <- readRDS(paste(data.loc,'meta_data.RDS',sep=''))
}else {
  source(paste(help.fun.loc,'/get_meta_data.R',sep=''))
  saveRDS(meta.data.neus,paste(data.loc,'meta_data.RDS',sep=''))
}
#------------------------------------------------------------------------------------------------------------------
#this is a mess, using two different codes to pull data
#------------------------------------------------------------------------------------------------------------------
all.new.england <- unlist(new_england_fossil)
site.index <- grep('dataset.id',names(unlist(new_england_fossil)))# this name is in fact wrong but I correctly call dataset.id 
site.name.index <- grep('site.name',names(unlist(new_england_fossil))) 
#the abbreviatiosn lat and lon occur twice, we therefore need some additional work to extract real lats and lons only
site.lat.index <- grep('lat',names(unlist(new_england_fossil))) 
site.lon.index <- grep('lon',names(unlist(new_england_fossil)))  
site.lat.acc.index <- grep('lat.acc',names(unlist(new_england_fossil))) 
site.lon.acc.index <- grep('long.acc',names(unlist(new_england_fossil))) 
site.lat.index.eff <- site.lat.index[!site.lat.index%in%site.lat.acc.index] #do not use the indexes that intersect
site.lon.index.eff <- site.lon.index[!site.lon.index%in%site.lon.acc.index]



#extract some meta data
site_ids <- as.numeric(all.new.england[site.index]) 
site.name <- all.new.england[site.name.index]
site.lat <- round(as.numeric(all.new.england[site.lat.index.eff]),2)
site.lon <- round(as.numeric(all.new.england[site.lon.index.eff]),2)

#--------------------------------------------------------------------------------------------------------------------
#fuse the two calls for data, we only use data from the US
index_us <- which(site_ids%in%meta.data.neus$datasetID)
site_ids_us <- site_ids[index_us]
site_ids_us <- sort(site_ids_us)

#--------------------------------------------------------------------------------------------------------------------
# Download data with a loop
#
#--------------------------------------------------------------------------------------------------------------------
if('all_data_NE.RDS'%in%list.files(data.loc)){
  all.data.NE <- readRDS(paste(data.loc,'all_data_NE.RDS',sep=''))
}else {
all.data.NE <- lapply(site_ids_us,function(x) {# run over the site indexes found above
  # download data for a specified site
  data.download <- get_download(x) 
  # get depth of a sample
  depth <- data.download[[1]]$sample.meta$depth
  lon <- data.download[[1]]$dataset$site$long
  lat <- data.download[[1]]$dataset$site$lat 
  name.site <- data.download[[1]]$dataset$site$site.name
  site.id <- data.download[[1]]$dataset$site.data$site.id
  #find chronologies of this core
  chronologies <- data.download[[1]]$chronologies
  #some cores have more than one chronology stored, we have to check if any of these chronologies contains a date in the 
  #period we are interested in 
  chronology_contain <- sapply(1:length(chronologies),function(zz) {
    sum((chronologies[[zz]]$age>=younger.age)&(chronologies[[zz]]<=older.age)[,'age'],na.rm=TRUE)})
  #if there is more than one chronology containing required dates, we take the one stored first in neotoma
  chronology_use <- min(which(chronology_contain>0))
  
  if(is.finite(chronology_use)) {
  chronology <- chronologies[[chronology_use]]$age
  #load all counts in the first place
  all.counts <- data.download[[1]]$counts
  #find Lycopdoium, Eucalyptus and Microspheres
  ind.lyco <- grep('Lycopodium',colnames(all.counts))
  ind.euca <- grep('Eucalyptus',colnames(all.counts))
  ind.micro <- grep('Micros',colnames(all.counts))
  ind.tot <- c(ind.lyco,ind.euca,ind.micro)
  if(length(ind.tot)>0) all.counts <- all.counts[,c(-ind.tot)]
  #load counts that fall withing a given bracket
  counts <- as.matrix(all.counts[which((chronology>=younger.age)&(chronology<=older.age)==TRUE),]) #Remove lycopodium Eucalyptus asf.
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
  depth.samples <- depth[which((chronology>=younger.age)&(chronology<=older.age)==TRUE)]
  # load chornology of selected data 
  chronology <- chronology[which((chronology>younger.age)&(chronology<older.age)==TRUE)] #for some weird reason I have to set this TRUE,
  #else it also takes depths with NA
   #return data for later use
  list(counts = counts, percentages = percentages,percentages.clean = percentages.clean,chronology = chronology,
       dataset_idx = x,site_idx = site.id,depth = depth.samples,lat = lat, lon=lon,site.name = name.site)
  }
})
saveRDS(all.data.NE,paste(data.loc,'all_data_NE.RDS',sep=''))
}

#----------------------------------------------------------------------------------------------------------------------




#site names

#find data that has no entry (i.e. no data in the requested period)
null_idx <- sapply(all.data.NE,function(x) !is.null(x))
#only retain sites that have entries
all.data.NE.cleaned <- lapply((1:length(all.data.NE))[null_idx],function(x) all.data.NE[[x]])


#create some meta data
dataset_id <- sapply(1:length(all.data.NE.cleaned),function(x) all.data.NE.cleaned[[x]]$dataset_idx)
site_id <- sapply(1:length(all.data.NE.cleaned),function(x) all.data.NE.cleaned[[x]]$site_idx)
lon.site <- sapply(1:length(all.data.NE.cleaned),function(x) all.data.NE.cleaned[[x]]$lon)
lat.site <- sapply(1:length(all.data.NE.cleaned),function(x) all.data.NE.cleaned[[x]]$lat)
name.site <- sapply(1:length(all.data.NE.cleaned),function(x) all.data.NE.cleaned[[x]]$site.name)
write.csv(data.frame(site.names= name.site,dataset.id = dataset_id,site.id = site_id,lon=lon.site,lat= lat.site),row.names=FALSE,paste0(wd,'expert_elicitation/data/site_names.csv'))




#-------------------------------------------------------------------------------------------------------------------
#species included for elicitation plots
#-------------------------------------------------------------------------------------------------------------------
spec_inc <- c('Tsuga','Picea','Pinus','Abies','Betula','Quercus','Fagus','Fraxinus','Castanea','Carpinus',
              'Cyperaceae','Poaceae','Artemisia','Asteraceae un','Rumex','Ambrosia')
              
#-------------------------------------------------------------------------------------------------------------------
# transform data into percentages
#------------------------------------------------------------------------------------------------------------------
percentages_spec <- sapply(1:length(all.data.NE.cleaned),function(y){#
  #extract the count data
  counts <- as.matrix(all.data.NE.cleaned[[y]]$counts) 
  #there are different formats for naming the pollen data  the if clauses ensure that we have properly named pollen
  if (length(colnames(counts))>0) column_names <- colnames(counts)
  if (length(colnames(counts))==0) column_names <- dimnames(counts)[[2]]
  # select specific taxa (note these are not cleaned, i.e. Pinus can have more than one entry)
  col_use <- unlist(sapply(spec_inc,function(x) grep(x,column_names))) 
  #select the specific taxa from the table
  counts_use <- as.matrix(counts[,col_use])
  #if as.matrix command transposes our data, we have to backtranspose it
  if(ncol(counts_use) == 1) perc_use <- t(counts_use)/sum(counts_use)
  #if there is more than one sample, add the sums in the samples
  if(min(dim(counts_use))>1) perc_use <- counts_use/rowSums(counts_use)
  #ensure a value is returned
  round(100*perc_use) 
})

#-------------------------------------------------------------------------------------------------------
fake.counts <- matrix(nrow=1,rep(0,length(spec_inc)))
colnames(fake.counts) <- spec_inc
fake.counts = as.data.frame(fake.counts)
#-------------------------------------------------------------------------------------------------------
# bring pollen to genus level
#-------------------------------------------------------------------------------------------------------
cleaned_percentages_spec <- lapply(1:length(all.data.NE.cleaned), function(y){ #length(percentages_spec),
   single.site <- sapply(spec_inc, function(x) {
    #for(i in 1:length(spec_inc)) { 
    # x <- spec_inc[i]
    #find elements that match one genus name
    pollen_perc <- as.data.frame(all.data.NE.cleaned[[y]]$percentages)
    pollen_perc <- replace(pollen_perc,is.na(pollen_perc),0)
   # if(is.numeric(nrow(pollen_perc))) pollen_perc <- pollen_perc[rowSums(pollen_perc)>0,]
    pollen_perc <-  Merge(pollen_perc,fake.counts,split=TRUE)$pollen_perc
    idx <- grep(x,colnames(pollen_perc))
    #if(is.null(colnames(counts_spec[[y]]))) idx <- grep(x,names(counts_spec[[y]]))
    if((nrow(pollen_perc)==1) & (length(idx) == 1)) cleaned_perc <- pollen_perc[idx]
    if((nrow(pollen_perc)==1) & (length(idx) > 1)) cleaned_perc <- sum(pollen_perc[idx])
    if((nrow(pollen_perc)>1) & (length(idx) == 1)) cleaned_perc <- pollen_perc[,idx]
    if((nrow(pollen_perc)>1) & (length(idx) > 1)) cleaned_perc <- rowSums(pollen_perc[,idx])
    cleaned_perc
  })
  if(is.list(single.site)) single.site <- matrix(ncol = length(spec_inc),as.numeric(single.site),byrow=TRUE)
  single.site
})







#-----------------------------------------------------------------------------------------------------------------------
#now we have to find the site names asf
#---------------------------------------------------------------------------------------------------------------------
site_ids_used <- site_ids[null_idx]
site_lon <- site.lon[null_idx]
site_lat <- site.lat[null_idx]
site_name <- site.name[null_idx]



#load modified strat plot
source(paste(help.fun.loc,'strat_plot_mod.R',sep=''))








pdf(paste(plot.loc,"/pollen_diagrams_NEUS_sqrt.pdf",sep=''),width = 11.69,height = 8.29)
  map(database = 'state',xlim=c(-85,-65),ylim=c(38,49))
  points(lon.site,lat.site,pch = 16,cex = 0.75)
  for (i in 1:length(cleaned_percentages_spec) ){
    layout(1)
    par(oma=c(4,4,2,2))
    if(is.numeric(nrow(cleaned_percentages_spec[[i]]))) {
      col_bar <- rep('grey',nrow(cleaned_percentages_spec[[i]]))
      if(nrow(cleaned_percentages_spec[[i]])>4) col_bar[seq(5,nrow(cleaned_percentages_spec[[i]]),5)] <- 'blue' 
      }
    else col_bar <- rep('grey',1)
  
  
    spec_plot <- apply(cleaned_percentages_spec[[i]],2,function(x) max(x)>2) 
    spec_plot[(length(spec_plot)-1):length(spec_plot)] <- rep(TRUE,2)
    # 2
    if(nrow(cleaned_percentages_spec[[i]])>1) strat.plot.mod(d =sqrt(cleaned_percentages_spec[[i]][,spec_plot]),
                                                                 yvar =   all.data.NE.cleaned[[i]]$depth,
                                                                 scale.percent=TRUE,col.bar = col_bar, 
                                                                 y.rev = TRUE,ylabel = 'Depth [cm]',
                                                                 title = paste(i,' ', name.site[i],'     Coord: ',
                                                                              round(abs(lon.site[i]),2), 'W ', 
                                                                              round(lat.site[i],2),'N',
                                                                              sep = ''),cex.title = 1.2 )

  }
dev.off()




#---------------------------------------------------------------------------------------------------------------------
#Store elicitation metadata

#find plotted sites 
ind.plot <- sapply(1:length(cleaned_percentages_spec),function(i){
  nrow(cleaned_percentages_spec[[i]])>1
})

data.for.csv <- data.frame(name.site[ind.plot],round(lon.site[ind.plot],2),round(lat.site[ind.plot],2))
colnames(data.for.csv) <- c('Site','Longitude','Latitude')
write.csv(data.for.csv,paste(data.loc,'/pollen_meta.csv',sep=''))
write.table(data.for.csv,paste(data.loc,'/site_information.txt',sep=''))
write.table(site_ids_us,paste(data.loc,'/site_ids.txt',sep=''))
write.table(which(ind.plot==TRUE),paste(data.loc,'/plot_index.txt',sep=''))

