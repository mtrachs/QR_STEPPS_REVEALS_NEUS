#-------------------------------------------------------------------------------------------------------------------------
#load data from elicitation exercise
#-------------------------------------------------------------------------------------------------------------------------
fil_names <- list.files('~/workflow_stepps_calibration/expert_elicitation/data/Pre-settlement Pollen/')
index_csv <- grep('.csv',fil_names)

#find sample indicated by experts
pre_settlement_sample <- 
sapply(index_csv,function(x) {
  file <- read.csv(paste('~/workflow_stepps_calibration/expert_elicitation/data/Pre-settlement Pollen/',fil_names[x],sep=''))
  as.numeric(as.character(file[,grep('settlement',colnames(file))]))
})

# find the number of agreements among the three scientist we have so far
# number_of_agreements <- apply(pre_settlement_sample,1,function(x) max(table(x)))
# 
# # find the sample number on which they agree
# sample_number <- apply(pre_settlement_sample,1,function(x)  ifelse(length(table(x))>0,
#   as.numeric(names(which.max(table(x)))),NA))
# 
# 
# #find lakes with agreement
# sites.use <- number_of_agreements>1
# 
# #find names of sites used for elicitation plots
# first.csv <- read.csv(paste('~/workflow_stepps_calibration/expert_elicitation/data/Pre-settlement Pollen/',
#                             fil_names[min(index_csv)],sep=''))
# names.elicitation <- first.csv[,'Site']#[,grep('Site',colnames(first.csv))]
# names.sites.use <- names.elicitation[sites.use] 
# 
# 
# 
# 
# 
# #combine sample number and site number
# site.sample <- cbind(1:123,sample_number)[sites.use,]
# colnames(site.sample)[1] <- 'site.index'
# site.sample <- as.data.frame(site.sample)
# 
# #---------------------------------------------------------------------------------------------------------------------
# #produce a random sample from the answers:
# 
# sample_use <-   replicate(100,{
#     apply(pre_settlement_sample,1,function(x) sample(x,1))
# })


#--------------------------------------------------------------------------------------------------------------------
#look at certainty of the analysts

sample_certainty <- 
  sapply(index_csv,function(x) {
    file <- read.csv(paste('~/workflow_stepps_calibration/expert_elicitation/data/Pre-settlement Pollen/',fil_names[x],sep=''))
    as.numeric(as.character(file[,grep('certainty',colnames(file))]))
  })


sample_certainty1 <-  apply(sample_certainty,1,function(x) sum(x <2,na.rm=TRUE))
table(sample_certainty1)



#-----------------------------------------------------------------------------------------------------------------------
#find names of sites 
site.information <- read.table('~/workflow_stepps_calibration/expert_elicitation/data/site_information.txt')

pre_settlement_sample1 <- t(apply(pre_settlement_sample,1,function(x) sort(x,na.last=TRUE)))

pre_settlement_sample1 <- cbind(site.information,pre_settlement_sample1)


# pre_settlement_sample1[sample_certainty1==0,]
# pre_settlement_sample1[sample_certainty1==1,]
# pre_settlement_sample1[sample_certainty1==2,]
# pre_settlement_sample1[sample_certainty1==3,]
# pre_settlement_sample1[sample_certainty1==4,] # all of these samples are currently in.. obviously not what we want


#check how many experts agree
# lapply(0:4,function(x) apply(pre_settlement_sample1[sample_certainty1==x,],1,function(x) max(table(x),na.rm=TRUE)))

#try to use all that have zero or one uncertain expert 
sites.certain <- lapply(0:1,function(x) pre_settlement_sample1[sample_certainty1==x,])
sites.index <- unlist(sapply(1:length(sites.certain),function(x) as.numeric(rownames(sites.certain[[x]]))))
samples.certain <- unlist(sapply(1:length(sites.certain),function(x) 
  sites.certain[[x]][,colnames(sites.certain[[x]])=='3']))

site.sample <- as.data.frame(cbind(sites.index,samples.certain))
colnames(site.sample) <- c("site.index","sample_number")
#looks strange because I am merging two groups, those with agreement of all experts and with one expert being uncertain 
#about the existence of a pre-settlement sample
site.sample <- site.sample[order(site.sample$site.index),]

