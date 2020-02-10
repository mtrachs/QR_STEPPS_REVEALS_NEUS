#--------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------
#devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)
library(dplyr)
library(DT)
library(neotoma)
library(sp)
library(fields)
library(rgdal)
library(abind)
library(rstan)

#-----------------------------------------------------------------------------------------------
setwd(paste0(wd,'prediction/'))
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'


load(paste0(wd,'vegetation/data/township_data_13_taxa_6796_cells_120_knots.rdata'))# do I need to load that file? Yes, gives some kind of spatial structure
spat.param.umw <- read.csv(paste0(wd,'vegetation/data/spatial_pars_UMW.csv'))
taxa <- c('Ash','Beech','Birch','Chestnut','Hemlock','Hickory','Maple',"Oak",'Other conifer',
          'Other hardwood','Pine','Spruce','Tamarack') 
taxa.umw <- c('Ash','Beech','Birch','Elm','Hemlock','Maple',"Oak",'Other conifer', 'Other hardwood','Pine','Spruce','Tamarack') 
spat.param.umw$taxon <- taxa.umw

rho <- spat.param.umw$rho[taxa.umw%in%taxa]
rho.neus <- rep(NA,length(taxa))
names(rho.neus) <- taxa
rho.neus[taxa%in%taxa.umw] <- rho
rho.neus[is.na(rho.neus)] <- median(rho.neus[!is.na(rho.neus)])
rho <- rho.neus


eta <- spat.param.umw$eta[taxa.umw%in%taxa]
eta.neus <- rep(NA,length(taxa))
names(eta.neus) <- taxa
eta.neus[taxa%in%taxa.umw] <- eta
eta.neus[is.na(eta.neus)] <- median(eta.neus[!is.na(eta.neus)])
eta <- eta.neus


#-----------------------------------------------------------------------------------------------

num_sites <- 110

load(paste0(wd,'calibration/data/for_prediction_110_sites_only_abies_new_species.RData'))

saveRDS(pollen_coords,paste(data.loc,'final_pollen_coords_calib.RDS',sep=''))
saveRDS(y,paste(data.loc,'settlement_era_pollen_final.RDS',sep=''))


N <- N_cells

#------------------------------------------------------------------------------------------------
#produce knots
#------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# have to load parameters from calibration model
#--------------------------------------------------------------------------------------------
library(rstan)
source(paste(help.fun.loc,'pred_helper_funs.r',sep='')) #make sure I use medians of a
source(paste(help.fun.loc,'process_funs.r',sep='')) #make sure I use medians of a

# in a first sep I will load phi and gamma from power law kernel Ka_Kgamma
source('R/build_cal_main.r') # this is strange....
runs1 <- list(runs[[4]])
run <- runs1[[1]]
kernel <- run$kernel
num_a <- run$one_a
one_psi <- run$one_psi
handle <- strsplit(run$suff_fit,'_A')[[1]][1]

fname <- paste0(wd,'results/data/cal_pl_Ka_Kgamma_EPs_modified_a_110.csv') 
fit <- read_stan_csv(fname)
post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)

param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('gamma','phi') #only take eta and rho
est.e.bayes <- colMeans(post[,1,param.e.bayes])
gamma.all <- est.e.bayes[grep('gamma',names(est.e.bayes))]
gamma <- gamma.all
phi.all <- est.e.bayes[grep('phi',names(est.e.bayes))]
phi <- phi.all

#this is certainly not correct yet
w <- build_weight_matrix(post = post,d = t(d),idx_cores = idx_cores,
                         N = N_cells,N_cores =  N_cores,run = run)

#w <- w/sum(w)# should not do that!

sum_w_pot <- build_sumw_pot(post = post, K=K,N_pot = N_pot,
                            d_pot =  d_pot,run = run)

if(length(gamma)==1) {
  gamma <- rep(gamma,K)
 w1 <- array(dim=c(K,N_cores,N_cells))
 for(i in 1:K) {w1[i,,]<-w}
 w <- w1
}
#-------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
res <- 1 # 
T <- 1
#-------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#change d_knots
colnames(knot_coords) <- c('east','north')
knot_coords <- as.data.frame(knot_coords)
veg_coords <- as.data.frame(veg_coords)
knot_coords <- knot_coords[((knot_coords$east>min(veg_coords$meters.east)) & (knot_coords$north>min(veg_coords$meters.north))),]

d_inter <- rdist(veg_coords,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6  
N_knots <- nrow(knot_coords)
#--------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
#store data in .dump file
stan_rdump(
  c('K', 'N','T', 'N_knots', 'N_cores',
    'y', 'res',
    'sum_w_pot',
    'rho','eta','gamma','phi',
    'idx_cores',
    'd','d_knots','d_inter','w'
    ), 
file=paste(paste0(wd,'vegetation/data/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'_sites_final_spat_param_umw.dump',sep="")))


save(K, N,T, N_knots, N_cores,
       y, res,
       sum_w_pot,
       rho,eta,gamma,phi,
       idx_cores,
       d,d_knots,d_inter,w,
      pollen_coords,r,veg_coords,
     #coords.agg.final
#, 
file=paste(paste0(wd,'vegetation/data/prediction_',K,'_taxa_',N_cells,'_cells_',N_knots,'_knots_',handle,'_',num_sites,'_sites_final_spat_param_umw.rdata',sep="")))


