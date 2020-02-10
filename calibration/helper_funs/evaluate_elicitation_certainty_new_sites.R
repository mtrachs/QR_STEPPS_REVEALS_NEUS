####################################################################################################################
#
####################################################################################################################
elicitation.loc <- paste0(wd,'/expert_elicitation/data_new/')

fil_names <- list.files(elicitation.loc)
index_csv <- grep('.csv',fil_names)

index.New.England <- grep('England',fil_names) 
index.NEUS <- grep('NEUS',fil_names) 


#find sample indicated by experts
pss.New_England <- 
  sapply(index.New.England,function(x) {
    file <- read.csv(paste(elicitation.loc,fil_names[x],sep=''))
    rn <- as.character(file[,'Site.name'])
    sample <- as.numeric(as.character(file[,grep('Sample',colnames(file))]))
    names(sample) <- rn
    if(x==min(index.New.England)) saveRDS(file[,c('lon','lat')],
                                          paste0(wd,'calibration/data/coordinates_New_England.RDS'))
    sample
  })

pss.NEUS <- 
  sapply(index.NEUS,function(x) {
    file <- read.csv(paste(elicitation.loc,fil_names[x],sep=''))
    rn <- as.character(file[,'X'])
    sample <- as.numeric(as.character(file[,grep('Sample',colnames(file))]))
    names(sample) <- rn
    if(x==min(index.NEUS)) saveRDS(file[,c('lon.site','lat.site')],
                                          paste0(wd,'calibration/data/coordinates_not_in_neotoma_NEUS.RDS'))
    sample
    })


