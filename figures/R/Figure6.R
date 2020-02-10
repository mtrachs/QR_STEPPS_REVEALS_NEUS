library(rstan)

setwd(paste0(wd,'figures/'))
output_location <- 'data/'
files_stan <- list.files(output_location)

summary_ppe <- lapply(1:length(files_stan), function(x) matrix(ncol=12, nrow = 3))

for(i in 1:length(files_stan)){
  posterior <- read_stan_csv(paste(output_location,files_stan[i],sep=''))
  samples <- posterior@sim$samples[[1]]
  index_use <- grep('phi',names(samples))
  ppe <- samples[-c(1:250),index_use]
  summary_ppe[[i]] <- ppe#round(apply(ppe,2,function(x) quantile(x,probs=c(0.025,0.5,0.975))))
  colnames(summary_ppe[[i]]) <- sort(c('Hemlock','Tamarack','Pine','Birch','Oak','Beech','Other Conifer','Spruce','Maple',
                                       'Ash','Elm','Other Hardwood'))
}

run_names <- matrix(ncol=2,unlist(strsplit(files_stan,'_A')),byrow=TRUE)[,1]

names(summary_ppe) <- run_names


summary_ppe_umw <- summary_ppe




#-------------------------------------------------------------------------------------------------------------------
#read productivity estimates for NEUS
#-------------------------------------------------------------------------------------------------------------------
#setwd('~/workflow_stepps_calibration/comparison_UMW_NEUS/')
source(paste0(wd,'comparison_UMW_NEUS/R/load_NEUS.r'))
setwd(paste0(wd,'figures/'))
figures_location <-'figures/'
summary_ppe_NEUS <- summary_ppe
names(summary_ppe_NEUS) <- runs_unlist1

#-----------------------------------------------------------------------------------------------------
n_sample <- 10000

ppe_difference <-
  lapply(runs_unlist1,function(y){
    umw <- summary_ppe_umw[[min(which(names(summary_ppe_umw)%in%y))]]
    neus <- summary_ppe_NEUS[[min(which(names(summary_ppe_NEUS)%in%y))]]
    neus <- neus[,colnames(neus)%in%colnames(umw)]
    taxa <- colnames(neus)
    
    ppe_difference <- 
      sapply(taxa,function(x){
        neus_sample <- sample(neus[,x],n_sample,replace = TRUE)
        umw_sample <- sample(umw[,x],n_sample,replace = TRUE)
        neus_sample - umw_sample
      })
  })

names(ppe_difference) <- runs_unlist1

differences_quantiles <-
  lapply(runs_unlist1,function(y){
    apply(ppe_difference[[y]],2,function(x) round(quantile(x,probs=c(0.025,0.5,0.975)),2))
  })

names(differences_quantiles) <- runs_unlist1


#---------------------------------------------------------------------------------------------------------------------
#Figure comparing productivity estimates
#---------------------------------------------------------------------------------------------------------------------
setwd(paste0(wd,'figures/'))
figures_location <-'figures/'
#take Ka_Kgamma and make a figure of estimates for both areas
taxa <- colnames(summary_ppe_umw[[1]])

run_names <- unique(names(summary_ppe_umw))
run_names1 <-unlist(strsplit(run_names,'cal_'))[seq(2,2*length(run_names),2)] 
run_names1 <-unlist(strsplit(run_names1,'_EPs'))

taxa.compare <- sort(c('Hemlock','Tamarack','Pine','Birch','Oak','Beech','Other Conifer','Spruce','Maple',
                       'Ash','Other Hardwood')) 

pdf(paste(figures_location,'Figure_6.pdf',sep=''),height=8,width = 10)
#sapply(run_names,function(rn){
rn <- "cal_pl_Ka_Kgamma_EPs" 

  ppe_umw <- apply(summary_ppe_umw[[min(which(names(summary_ppe_umw)==rn))]],2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
  ppe_umw <- as.data.frame(t(ppe_umw))
  ppe_umw <- as.data.frame(ppe_umw[rownames(ppe_umw)%in%taxa.compare,])
  ppe_neus <- apply(summary_ppe_NEUS[[which(names(summary_ppe_NEUS)==rn)]],2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
  ppe_neus <- as.data.frame(t(ppe_neus))
  ppe_neus <- ppe_neus[rownames(ppe_neus)%in%taxa.compare,]
  
  par_stats_umw = data.frame(name=taxa.compare, mu=ppe_umw$`50%`, lb=ppe_umw$`2.5%`, ub=ppe_umw$`97.5%`, handle=rep('UMW', length(taxa.compare)))
  par_stats_neus = data.frame(name=taxa.compare, mu=ppe_neus$`50%`, lb=ppe_neus$`2.5%`, ub=ppe_neus$`97.5%`, handle=rep('NEUS', length(taxa.compare)))
  
  dat = rbind(par_stats_umw,par_stats_neus)
  colnames(dat)[5] <- 'Region'
  
  p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=Region, colour=Region)) + 
    geom_point(size=4, position=position_dodge(width=0.5)) + 
    geom_errorbar(aes(ymin=lb, ymax=ub), width=.2, position=position_dodge(width=0.5)) + 
    xlab("Taxon") + ylab(parse(text='phi')) +
    coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                      axis.title.y=element_text(size=20), 
                                      axis.text.x=element_text(size=rel(1.3)),
                                      axis.text.y=element_text(size=rel(1.3)))
  
  print(p)
  #print(p + ggtitle(run_names1[run_names==rn]))
  
#})
dev.off()

