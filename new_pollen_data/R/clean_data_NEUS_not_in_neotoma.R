################################################################################################################################
# load additional pollen data not in neotoma
################################################################################################################################
setwd(paste0(wd,'new_pollen_data/'))
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

pollen.names.cleaned <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  cn <- colnames(pollen.counts)
  cn <- gsub('J.', 'Juglans',cn,fixed=TRUE)
  cn <- gsub('A.', 'Alnus',cn,fixed=TRUE)
  cn <- gsub('F.', 'Fraxinus',cn,fixed=TRUE)
  ind <- grep('P.',cn,ignore.case=FALSE,fixed = TRUE)
  if(length(ind)<1){
    cn}else{
      Pinus.index <- min(grep('Pinus',cn))
      Picea.index <- min(grep('Picea',cn))
      taxon.first <- which.min(c(Pinus.index,Picea.index))
      pinaceae.index <- sort(c(Pinus.index,Picea.index))
      
      ind1 <- ind[((ind>pinaceae.index[1])&(ind<pinaceae.index[2]))]
      ind2 <- ind[ind>pinaceae.index[2]]
      if(taxon.first==1){
        cn[ind1] <- 'Pinus'
        cn[ind2] <- 'Picea'
      }else{
        cn[ind1] <- 'Picea'
        cn[ind2] <- 'Pinus' 
      }
    }
  cn
})


pollen.names.abb <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  grep('J.',colnames(pollen.counts),fixed=TRUE)
})

pollen.counts.cleaned <- sapply(files[pollen.files],function(x){
  pollen.counts <- read.csv(paste(data.loc,x,sep=''))
  cn <- colnames(pollen.counts)
  #clean colnames
  cn <- gsub('J.', 'Juglans',cn,fixed=TRUE)
  cn <- gsub('A.', 'Alnus',cn,fixed=TRUE)
  cn <- gsub('F.', 'Fraxinus',cn,fixed=TRUE)
  ind <- grep('P.',cn,ignore.case=FALSE,fixed = TRUE)
  if(length(ind)<1){
    cn}else{
      Pinus.index <- min(grep('Pinus',cn))
      Picea.index <- min(grep('Picea',cn))
      taxon.first <- which.min(c(Pinus.index,Picea.index))
      pinaceae.index <- sort(c(Pinus.index,Picea.index))
      
      ind1 <- ind[((ind>pinaceae.index[1])&(ind<pinaceae.index[2]))]
      ind2 <- ind[ind>pinaceae.index[2]]
      if(taxon.first==1){
        cn[ind1] <- 'Pinus'
        cn[ind2] <- 'Picea'
      }else{
        cn[ind1] <- 'Picea'
        cn[ind2] <- 'Pinus' 
      }
    }
  colnames(pollen.counts) <- cn
  
  ind.halfs <- grep('half',cn)  
  if(length(ind.halfs)>0) pollen.counts[,ind.halfs] <- pollen.counts[,ind.halfs]/2
  ind.tot <- grep('tot',cn,ignore.case=TRUE)
  if(length(ind.tot)>0){
    pollen.counts <- pollen.counts[,-c(ind.tot)]
    cn <- colnames(pollen.counts)
  } 
  ind.sum <- grep('sum',cn,ignore.case=TRUE)
  if(length(ind.sum)>0){
    pollen.counts <- pollen.counts[,-c(ind.sum)]
    cn <- colnames(pollen.counts)
  }
  ind.stomata <- grep('stom',cn,ignore.case=TRUE)
  if(length(ind.stomata)>0){
    pollen.counts <- pollen.counts[,-c(ind.stomata)]
    cn <- colnames(pollen.counts)
  }
  ind.depth <- grep('depth',cn,ignore.case=TRUE)
  if(length(ind.depth)>0){
    pollen.counts <- pollen.counts[,-c(ind.depth)]
    cn <- colnames(pollen.counts)
  }
  ind.sample <- grep('Sample',cn,ignore.case=TRUE)
  if(length(ind.sample)>0){
    pollen.counts <- pollen.counts[,-c(ind.sample)]
    cn <- colnames(pollen.counts)
  }
  ind.core <- grep('core',cn,ignore.case=TRUE)
  if(length(ind.core)>0){
    pollen.counts <- pollen.counts[,-c(ind.core)]
    cn <- colnames(pollen.counts)
  }
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
site_names <- matrix(unlist(strsplit(site_names,'_')),ncol=2,byrow=TRUE)[,1]
saveRDS(site_names,paste(data.loc,'additional_pollen_site_names.RDS',sep=''))

