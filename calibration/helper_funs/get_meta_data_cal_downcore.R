library(maps)

#get_neo <- function(){
 
  
  gpids <- get_table(table.name='GeoPoliticalUnits')

  states = c('Maine','New Hampshire','Vermont','New York','Massachusetts','Rhode Island',
             'Connecticut','Pennsylvania','New Jersey')
  
  gp_rows <- vector(length = length(states))
  
  for (x in 1:length(states)) {
    gp_rows[x] <-  which(gpids$GeoPoliticalName == states[x])
  }
  
  
  
  ID <- vector(length=length(states))
  
  for(x in 1:length(states)) {
    ID[x] <- gpids[which(gpids$GeoPoliticalName == states[x]),1]
  }
  
  
  meta <- get_dataset(datasettype='pollen', gpid=ID, ageyoung = -50,ageold=15000)


  meta
#}

#get_neo()



# get_calcote <- function(){
#   
#   clh.sites  <- read.csv('data/hotchkiss_lynch_calcote_meta.csv', stringsAsFactors = FALSE)
#   clh.counts <- read.csv('data/hotchkiss_lynch_calcote_counts.csv', stringsAsFactors = FALSE)
#   
#   clh.sites$name <- gsub(" ","", clh.sites$name, fixed=TRUE)
#   clh.counts$name <- gsub(" ","", clh.counts$name, fixed=TRUE)
#   
#   n <- nrow(clh.sites)
#   
#   ids    <- vector(mode="numeric", length=0)
#   states <- vector(mode="numeric", length=0)
#   lat    <- vector(mode="numeric", length=0)
#   long   <- vector(mode="numeric", length=0)
#   site   <- vector(mode="character", length=0)
#   PI     <- vector(mode="character", length=0)
#   descriptor <- vector(mode="character", length=0)
#   
#   site.count = 0
#   
#   for (i in 1:n){
#     
#     site.i <- as.character(clh.sites$name[i])
#     idx    <- which(clh.counts$name == site.i)
#     type   <- clh.counts[idx,1]
#     
#     if (length(type) > 2){
#       
#       site.count = site.count + 1
#       
#       ids    = c(ids, paste('CLH', site.count, sep=''))
#       states = c(states, 'wisconsin')
#       site   = c(site, site.i)
#       lat    = c(lat, clh.sites$lat[i])
#       long   = c(long, clh.sites$long[i])
#       PI     = c(PI, clh.counts$analyst[idx[1]])
#       
#     }
#     
#   }
#   
#   empty = rep(NA, length(ids))
#   
#   pollen_meta <- data.frame(datasetID=ids, handle=toupper(site), site=site, long=long, lat=lat, state=states, 
#                             pi=PI, pre=empty, settlement=empty, notes=empty, 
#                             stringsAsFactors=FALSE)
#   
#   # combine the meta and counts into a single table
#   pollen_all = NULL
#   for (i in 1:nrow(pollen_meta)){
#     idx <- which(clh.counts$name == pollen_meta$site[i])
#     
#     counts.i <- clh.counts[idx,!(colnames(clh.counts) %in% c('name', 'analyst'))]
#     
#     meta.i <- pollen_meta[rep(i, length(idx)),]
#     
#     pollen_all = rbind(pollen_all, data.frame(meta.i, counts.i, row.names=NULL))
#   }
#   
#   write.table(pollen_all, 'data/pollen_counts_calcote.csv', row.names=FALSE, sep=',')
#   
#   return(pollen_meta)
#   
# }
# 
# pollen_meta_neo     <- get_datasetids()
# pollen_meta_calcote <- get_calcote()
# 
# pollen_meta <- rbind(pollen_meta_neo, pollen_meta_calcote)
# pollen_meta <- pollen_meta[with(pollen_meta, order(handle)),]
# 
# write.table(pollen_meta, paste('pollen_meta_', Sys.Date(), '.csv', sep=''), row.names=FALSE, sep=',', na='')
