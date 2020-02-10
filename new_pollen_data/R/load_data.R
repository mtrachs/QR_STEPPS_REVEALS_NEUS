############################################################################################################################
# read additional data by Paciorek and McLachlan
############################################################################################################################
library(readr)
library(dplyr)
library(analogue)
library(stepps)
###########################################################################################

setwd(wd)
data.loc <- 'new_pollen_data/'
help.fun.loc <- 'calibration/helper_funs/' 
save.loc <- 'calibration/data/'
#############################################################################################################################
#load data of 11 sites 
# for all these data pinus total, picea total asf was set to 0
############################################################################################################################
pollen.data <- read.csv(file = paste(data.loc, 'Paciorek_McLachlan_additional.csv',sep=''))
lake.names <- as.character(pollen.data[seq(1,nrow(pollen.data),3),1]) 
taxa.names <- pollen.data[seq(2,nrow(pollen.data),3),]
taxon.names <- as.character(unique(unlist(taxa.names)))
pollen.counts <- pollen.data[seq(3,nrow(pollen.data),3),]

l.name <- 'Paciorek_McLachlan'

#make taxa translation table 
write.csv(sort(taxon.names),paste(data.loc,l.name,'taxon_names.csv',sep=''))
source(paste(help.fun.loc,'taxa_translation_only_abies_lake_names.R',sep=''))
pol_table <- readr::read_csv(paste(data.loc,l.name,'taxon_names_translated.csv',sep=''))


#########################################################################################################################
fake.counts <- matrix(nrow=1,rep(0,length(sort(unique(taxa.english)))))
colnames(fake.counts) <- sort(unique(taxa.english))
fake.counts <- as.data.frame(fake.counts)


calib_trans <- 
  sapply(1:length(lake.names),function(zz){
        daten <- t(as.matrix(unlist(pollen.data[seq(3,nrow(pollen.data),3)[zz],])))
        daten <- as.numeric(daten)
        daten[is.na(daten)] <- 0
        daten <- matrix(nrow=1,daten)
        colnames(daten) <- as.character(unlist(pollen.data[seq(2,nrow(pollen.data),3)[zz],]))
        daten <- as.data.frame(daten)
        agg.daten <- translate_taxa(daten, 
                                    pol_table,
                                    id_cols = colnames(daten)[1:2])

        joined.data <- analogue::join(agg.daten,fake.counts)$agg.daten
        pollen.counts <- joined.data[-c(1:2)]
        pollen.counts <- pollen.counts[order(colnames(pollen.counts))]
        joined.data <-c(joined.data[1:2],pollen.counts)
  })

additional_pollen <- t(calib_trans)


saveRDS(additional_pollen,paste(save.loc,'data_paciorek_mclachlan.RDS',sep=''))
