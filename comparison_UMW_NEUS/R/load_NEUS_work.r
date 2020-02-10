library(rstan)
library(ggplot2)
library(fields, quietly=TRUE)
library(abind)


#setwd(paste0(wd,'/'))
#wd = getwd()

#####################################################################################
# user pars
#####################################################################################

path_utils = paste0(wd,'results/utils/')
path_data  = paste0(wd,'results/data/')
path_out   = 'data/'#'output/'#'data/'
path_figs  = 'plots'

suff  =''
rerun = TRUE

kernel   = run$kernel
suff_fit = run$suff_fit
suff_dat = run$suff_dat

if (kernel == 'gaussian'){
  one_psi    = run$one_psi
  one_gamma  = run$one_gamma
  EPs        = run$EPs
} else if (kernel == 'pl'){
  one_a      = run$one_a
  one_b      = run$one_b
  one_gamma  = run$one_gamma
  EPs        = run$EPs
}

save_plots = TRUE
rescale    = 1e6

#####################################################################################
# read in data and source utils
#####################################################################################

source(file.path(path_utils, 'process_funs.r'))
source(file.path(path_utils, 'plot_funs.r'))
  
path_figs1 = sprintf('%s/%s', path_figs, suff_fit)
if (!file.exists(path_figs1)){
  dir.create(file.path(path_figs1))
}

#load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))
load(paste0(wd,'calibration/data/elicitation_neus_certainty_median_110_sites_only_abies_new_species.RData'))
taxa <- colnames(y)

file.name <- strsplit(suff_fit,'NEUS')
file.name <- file.name[[1]][1]
if(kernel=='pl'){
  file.name <- paste(file.name,'modified_a_',sep='')#ifelse(nchar(file.name)>6,paste('output_',file.name,sep=''),'output_neus_')
}





  fname = paste(path_data,file.name,'110.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  fit <- read_stan_csv(fname)
  post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)




samples <- post[,1,]
ppe <- samples[,grep('phi',colnames(samples))]
summary_ppe[[which(runs_unlist%in%suff_fit)]] <- ppe#round(apply(ppe,2,function(x) quantile(x,probs=c(0.025,0.5,0.975))))
colnames(summary_ppe[[which(runs_unlist%in%suff_fit)]]) <- sort(c('Hemlock','Tamarack','Pine','Birch','Oak','Beech','Other Conifer','Spruce','Maple',
                                     'Ash','Other Hardwood','Chestnut','Hickory'))


