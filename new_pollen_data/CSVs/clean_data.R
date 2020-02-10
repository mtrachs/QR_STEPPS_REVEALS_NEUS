################################################################################################################################
# load additional pollen data not in neotoma
################################################################################################################################
setwd('~/workflow_stepps_calibration/new_pollen_data/')
data.loc <- 'CSVs/'

files <- list.files(data.loc)
files <- files[grep('.csv',files)]

pollen.files <- grep('pol',files,ignore.case=TRUE)

#files[pollen.files]

#for(i in files[pollen.files]) {
pollen.names <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  colnames(pollen.counts)
})


pollen.names.abb <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  grep('J.',colnames(pollen.counts),fixed=TRUE)
})

pollen.counts.cleaned <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  cn <- colnames(pollen.counts)
  ind.halfs <- grep('half',cn)  
  if(length(ind.halfs)>0) pollen.counts[,ind.halfs] <- pollen.counts[,ind.halfs]/2
  ind.tot <- grep('tot',cn)
  if(length(ind.tot)>0) pollen.counts <- pollen.counts[,-c(ind.tot)]
  ind.sum <- grep('sum',cn)
  if(length(ind.sum)>0)pollen.counts <- pollen.counts[,-c(ind.sum)]
  ind.stomata <- grep('stomata',cn)
  if(length(ind.stomata)>0)pollen.counts <- pollen.counts[,-c(ind.stomata)]
  ind.depth <- grep('depth',cn,ignore.case=TRUE)
  if(length(ind.depth)>0)pollen.counts <- pollen.counts[,-c(ind.depth)]
  ind.sample <- grep('Sample',cn,ignore.case=TRUE)
  if(length(ind.sample)>0)pollen.counts <- pollen.counts[,-c(ind.sample)]
  ind.core <- grep('core',cn,ignore.case=TRUE)
  if(length(ind.core)>0)pollen.counts <- pollen.counts[,-c(ind.core)]
  pollen.counts <- replace(pollen.counts,is.na(pollen.counts),0)
  pollen.counts
})

saveRDS(pollen.counts.cleaned,paste(data.loc,'additional_pollen_clean.RDS',sep=''))


depths <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  cn <- colnames(pollen.counts)
  ind.depth <- grep('depth',cn,ignore.case=TRUE)
  ind.sample <- grep('Sample',cn,ignore.case=TRUE)
  if((length(ind.depth)>0)&(length(ind.sample)==0)) pollen.counts <- pollen.counts[,ind.depth]
  if((length(ind.depth)>0)&(length(ind.sample)>0)) pollen.counts <- pollen.counts[,ind.depth]
  if((length(ind.depth)==0)&(length(ind.sample)>0)) pollen.counts <- pollen.counts[,ind.sample]
  depth <- pollen.counts
  col.idx <- ncol(depth)
  if(length(col.idx)>0) depth <- depth[,2]
  if(length(grep('-',depth))>0){
    #depth <- gsub(' ','',depth)
    depth <- as.numeric(matrix(ncol=2,byrow=TRUE,unlist(strsplit(as.character(depth),'-')))[,1])
  }
  depth
})

saveRDS(depths,paste(data.loc,'additional_pollen_depths.RDS',sep=''))


site_names <- names(depths)
saveRDS(site_names,paste(data.loc,'additional_pollen_site_names.RDS',sep=''))

