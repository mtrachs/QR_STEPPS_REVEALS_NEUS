#######################################################################################################################
# load samples from expert elicitation
#######################################################################################################################

#sort samples
pss.NEUS  
pss.New_England
#######################################################################################################################
# Create index of sites for to be removed from NEUS list
# Benson and Doe have complete records in New England list
# From the Adirondacks, Jack wants to remove Queer, Windfall, Upper Wallface Pond, Lake Tear of the Clouds, 
# I want to remove Lake Arnold
# also remove lower part of Bloomingdale
#######################################################################################################################
removal.list <-  c('Benson','Doe','Queer','Windfall','Wallface', 'Tear','Arnold','lower')

remove.lakes <- 
  sapply(removal.list,function(x){
    grep(x,rownames(pss.NEUS))
  })

pss.NEUS1 <- as.matrix(pss.NEUS[-remove.lakes,])

coord.neus.not.neotoma <- readRDS(paste0(wd,'/calibration/data/coordinates_not_in_neotoma_NEUS.RDS'))
coord.neus.not.neotoma <- coord.neus.not.neotoma[-remove.lakes,]
saveRDS(coord.neus.not.neotoma,paste0(wd,'calibration/data/coordinates_not_in_neotoma_NEUS.RDS'))
#######################################################################################################################
#load pollen data associated with the new elicitation also load depths
#######################################################################################################################
pollen.NEUS_not_neotoma <- readRDS(paste0(wd,'new_pollen_data/CSVs/additional_pollen_clean.RDS'))
depths.NEUS_not_neotoma <- readRDS(paste0(wd,'new_pollen_data/CSVs/additional_pollen_depths.RDS'))
pollen.New.England <- readRDS(paste0(wd,'new_pollen_data/CSV_new/additional_pollen_New_England_clean.RDS'))
depths.New.England <- readRDS(paste0(wd,'new_pollen_data/CSV_new/additional_pollen_New_England_depths.RDS'))


pollen.NEUS_not_neotoma.clean <- sapply((1:length(pollen.NEUS_not_neotoma))[-remove.lakes],function(x){
  pollen.NEUS_not_neotoma[[x]]
})

names(pollen.NEUS_not_neotoma.clean) <- rownames(pss.NEUS1) 

depths.NEUS_not_neotoma.clean <- sapply((1:length(depths.NEUS_not_neotoma))[-remove.lakes],function(x){
  depths.NEUS_not_neotoma[[x]]
})

names(depths.NEUS_not_neotoma.clean) <- rownames(pss.NEUS1) 

########################################################################################################################
# choose the appropriate sample
########################################################################################################################
sample.use.NEUS <- apply(pss.NEUS1,1, function(x)
                      if(sum(is.na(x))>0){
                        max(x,na.rm=TRUE)
                      }else{
                        median(x)
                      })

settlement.era.sample.NEUS <- 
  lapply(1:length(sample.use.NEUS),function(x){
    sample.use <- sample.use.NEUS[x]      
    pollen.NEUS_not_neotoma.clean[[x]][sample.use,]
  })

names(settlement.era.sample.NEUS) <- rownames(pss.NEUS1) 

sample.use.New.England <- apply(pss.New_England,1, function(x)
  if(sum(is.na(x))>0){
    min(x,na.rm=TRUE) #take the minimum of guilder
  }else{
    median(x)
  })

settlement.era.sample.New.England <- 
  lapply((1:length(sample.use.New.England)),function(x){
    sample.use <- sample.use.New.England [x]      
    pollen.New.England[[x]][sample.use,]
  })

names(sample.use.New.England) <- rownames(pss.New_England) 
######################################################################################################################
settlement.era.depth.NEUS <- 
  sapply(1:length(sample.use.NEUS),function(x){
    sample.use <- sample.use.NEUS[x]      
    c(sample.use,depths.NEUS_not_neotoma.clean[[x]][sample.use])
  })

colnames(settlement.era.depth.NEUS) <- rownames(pss.NEUS1) 

settlement.era.depth.New.England <- 
  sapply((1:length(sample.use.New.England)),function(x){
    sample.use <- sample.use.New.England [x]      
    c(sample.use,depths.New.England[[x]][sample.use])
  })

colnames(settlement.era.depth.New.England) <- rownames(pss.New_England) 


write.csv(rbind(t(settlement.era.depth.NEUS),t(settlement.era.depth.New.England)),
          paste0(wd,'calibration/data/sample_nr_depth_new_data.csv'))



##########################################################################################################################
new_pollen_data <- c(settlement.era.sample.NEUS,settlement.era.sample.New.England)
names(new_pollen_data) <- c(names(settlement.era.sample.NEUS),names(sample.use.New.England))

saveRDS(new_pollen_data,paste0(wd,'calibration/data/settlement_era_data_not_in_neotoma.RDS'))
