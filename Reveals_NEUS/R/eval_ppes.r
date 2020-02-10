#######################################################################################################################
#standardize pollen prductivity estimates
#######################################################################################################################
library(reshape2)
library(ggplot2)

setwd(paste0(wd,'Reveals_NEUS'))
data.loc <- 'data/'

PPEs = read.csv(paste(data.loc,'PPEs.csv',sep=''))
PPEs = PPEs[which(!is.na(PPEs$ppe)),]
PPEs = PPEs[which(PPEs$use == 'y'),] 

taxon_list <- unique(PPEs$taxon)
# errors <- PPEs$error

PPEs$errors.estimated <- NA



#deal with missing errors, assume coefficient of variations is constant within taxon

for(tax in taxon_list){  
  estimate <- PPEs$ppe[PPEs$taxon==tax]  
  error <- PPEs$error[PPEs$taxon==tax] 
  error.prop <- error[!is.na(error)]/estimate[!is.na(error)]
  error.prop.use <- mean(error.prop[error.prop>0])
  error[is.na(error)] <- estimate[is.na(error)] *error.prop.use
  PPEs$errors.estimated[PPEs$taxon==tax] <- error
}

# unique(PPEs$taxon[!is.finite(PPEs$errors.estimated)])

# PPEs[PPEs$taxon=='Larch',]
# PPEs[PPEs$taxon=='Hickory',]
# PPEs[PPEs$taxon=='Ash',]

###########################################################################################################################3
# for taxa without any error estimate I use a global coefficient of variation 
error.proportions <- PPEs$error[!is.na(PPEs$error)]/PPEs$ppe[!is.na(PPEs$error)]
error.proportions.use <- mean(error.proportions[error.proportions>0])
PPEs$errors.estimated[!is.finite(PPEs$errors.estimated)] <-  PPEs$ppe[!is.finite(PPEs$errors.estimated)]* error.proportions.use 
PPEs1 <- PPEs
##################################################################################################################
#standardize with oak
dataset.id <- unique(PPEs$dataset)

dataset.id <- dataset.id[!dataset.id%in%c(11,13,28,29)]

PPEs$ppe.stand <- NA
PPEs$error.stand <- NA

for(i in dataset.id) {
  dat.use <- PPEs[PPEs$dataset==i,]
  PPEs$ppe.stand[PPEs$dataset==i] <- dat.use$ppe/dat.use$ppe[dat.use$taxon=='Oak']
  PPEs$error.stand[PPEs$dataset==i] <- dat.use$errors.estimated/dat.use$ppe[dat.use$taxon=='Oak']
  PPEs$error.stand[PPEs$dataset==i] <- sqrt(PPEs$ppe.stand[PPEs$dataset==i]^2*
                                              ((dat.use$errors.estimated[dat.use$taxon=='Oak']/dat.use$ppe[dat.use$taxon=='Oak'])^2+
                                                   (dat.use$errors.estimated/dat.use$ppe)^2))
}


PPEs <- PPEs[PPEs$dataset%in%dataset.id,]
write.csv(PPEs,paste(data.loc,'all_PPEs_considered.csv',sep=''))


PPEs_means <- aggregate(x=PPEs$ppe.stand,by = list(PPEs$taxon),FUN=mean)
PPEs_sds <- aggregate(x=PPEs$error.stand,by = list(PPEs$taxon),FUN=function(x) sqrt(mean(x[x>0]^2)))
PPEs_final <- cbind(PPEs_means$x,PPEs_sds$x)
rownames(PPEs_final) <- PPEs_means$Group.1
colnames(PPEs_final) <- c('ppe','ppe.error')

##################################################################################################################################
#if there are more than 2 values for taxon exclude the value with largest absolute distance to the mean
final_taxa <- sort(unique(PPEs$taxon))

PPE_Final <- data.frame(taxon = final_taxa,
                        ppe = rep(NA,length(final_taxa)),
                        ppe.error = rep(NA,length(final_taxa)),
                        ppe.sim = rep(NA,length(final_taxa)),
                        ppe.sim.error = rep(NA,length(final_taxa)))

for(i in final_taxa){
  taxon.index <- which(PPEs$taxon==i)
  
  if(length(taxon.index)==1) {
    PPE_Final$ppe[PPE_Final$taxon==i] <- PPEs$ppe.stand[taxon.index]
    PPE_Final$ppe.error[PPE_Final$taxon==i] <- PPEs$error.stand[taxon.index]
    simulation <-rnorm(10000,PPEs$ppe.stand[taxon.index],PPEs$error.stand[taxon.index])
    simulation <- simulation[simulation>0]
    PPE_Final$ppe.sim[PPE_Final$taxon==i] <- median(simulation)
    PPE_Final$ppe.sim.error[PPE_Final$taxon==i] <- sd(simulation)
  }
  if((length(taxon.index)==2 | (i=='Oak'))) {
    PPE_Final$ppe[PPE_Final$taxon==i] <- mean(PPEs$ppe.stand[taxon.index])
    PPE_Final$ppe.error[PPE_Final$taxon==i] <- sqrt(sum((PPEs$error.stand[taxon.index][PPEs$error.stand[taxon.index]>0])^2))
    simulation <- sapply(1:2,function(x) rnorm(10000,PPEs$ppe.stand[taxon.index][x],PPEs$error.stand[taxon.index][x]))
    simulation <- simulation[simulation>0]
    PPE_Final$ppe.sim[PPE_Final$taxon==i] <- median(simulation)
    PPE_Final$ppe.sim.error[PPE_Final$taxon==i] <- sd(simulation)
  }
  if(((length(taxon.index)>2) & (i!='Oak'))) {
    mittel <- mean(PPEs$ppe.stand[taxon.index])
    exclude <- which.max(abs(PPEs$ppe.stand[taxon.index]- mittel))
    taxon.index <- taxon.index[-exclude]
    PPE_Final$ppe[PPE_Final$taxon==i] <- mean(PPEs$ppe.stand[taxon.index])
    PPE_Final$ppe.error[PPE_Final$taxon==i] <- sqrt(sum((PPEs$error.stand[taxon.index][PPEs$error.stand[taxon.index]>0])^2))
    simulation <- sapply(1:length(taxon.index),function(x) rnorm(10000,PPEs$ppe.stand[taxon.index][x],PPEs$error.stand[taxon.index][x]))
    simulation <- simulation[simulation>0]
    PPE_Final$ppe.sim[PPE_Final$taxon==i] <- median(simulation)
    PPE_Final$ppe.sim.error[PPE_Final$taxon==i] <- sd(simulation)
  }
}


##################################################################################################################
#Alder only found in dataset 29 
# Andria standardized alder with mean of Birch
##################################################################################################################
alder_dataset <- PPEs1[PPEs1$dataset==29,]
write.csv(alder_dataset,paste(data.loc,'PPEs_alder_dataset.csv',sep=''))


alder.estimate_error <- alder_dataset$error[alder_dataset$taxon=='Alder']/alder_dataset$ppe[alder_dataset$taxon=='Alder']
birch.estimate_error <- alder_dataset$error[alder_dataset$taxon=='Birch']/alder_dataset$ppe[alder_dataset$taxon=='Birch']
birch.alder.ratio <- alder_dataset$ppe[alder_dataset$taxon=='Birch']/alder_dataset$ppe[alder_dataset$taxon=='Alder']  
uncertainty.stand.alder <- sqrt((1/birch.alder.ratio)^2*(alder.estimate_error^2 + birch.estimate_error^2)) 



alder.ppe.final <- PPE_Final$ppe[PPE_Final$taxon=='Birch']/birch.alder.ratio

alder.estimate <- data.frame(taxon = 'Alder',
                             ppe = alder.ppe.final,
                             ppe.error =  sqrt(alder.ppe.final^2*
                                                 ((uncertainty.stand.alder/alder.ppe.final)^2 + 
                                                 PPE_Final$ppe.sim.error[PPE_Final$taxon=="Oak"]^2))
                               ,
                             ppe.sim = alder.ppe.final,
                             ppe.sim.error = sqrt(alder.ppe.final^2*
                                                    ((uncertainty.stand.alder/alder.ppe.final)^2 + 
                                                       PPE_Final$ppe.sim.error[PPE_Final$taxon=="Oak"]^2)))

PPE_Final <- rbind(PPE_Final,alder.estimate)
PPE_Final <- PPE_Final[order(PPE_Final$taxon),]
colnames(PPE_Final) <- c('taxon','ppe.analytical','error.analytical','ppe','error')

saveRDS(PPE_Final,'data/PPEs_agg.RDS')





  

