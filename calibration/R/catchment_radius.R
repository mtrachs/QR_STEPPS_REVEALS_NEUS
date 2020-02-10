library(rstan)
library(abind)
library(reshape)

#

setwd(wd)
data.loc <- 'results/data/'
source('results/utils/capture_radius.r')
source('prediction/R/build_cal_main.r') # this is strange....
source('results/utils/process_funs.r')
load('calibration/data/elicitation_neus_certainty_median_110_sites_only_abies_new_species.RData')
taxa <- colnames(y)


dispersal_distance <- list(VAR_G=0,VAR_PL=0)

for (ii in 3:4){
    run <- runs[[ii]]
    kernel <- run$kernel
    num_a <- run$one_a
    one_psi <- run$one_psi
    run$handle <- strsplit(run$suff_fit,'_A')[[1]][1]

    #look at that again...
     # would have to change this if not run_pl
      if(kernel=='pl'){
       
          fname = paste0(wd,'results/data/cal_pl_Ka_Kgamma_EPs_modified_a_110.csv')
      }
      if(kernel=='gaussian'){
          fname = paste0('results/data/cal_g_Kpsi_Kgamma_EPs_110.csv')
      }
      fit <- read_stan_csv(fname)
      post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
    
    #taxa <- colnames(y)
 

radius = seq(8000,700000, by=4000)
    
    
disp_dat = dispersal_cdf(radius=radius, run=run, taxa=taxa, path_out=,d_pot=d_pot,post=post,rescale=10^6)

disp_mean = disp_dat$dat
disp_all  = disp_dat$r_all
#plot_dispersal_cdfs(disp_mean)


capture_25 = capture_radius(disp_mean, 0.25)
capture_50 = capture_radius(disp_mean, 0.5)
capture_70 = capture_radius(disp_mean, 0.7)
capture_90 = capture_radius(disp_mean, 0.9)

capture_50_ci = capture_radius_ci(disp_all, 0.5, radius,taxa= taxa)
capture_70_ci = capture_radius_ci(disp_all, 0.7, radius,taxa= taxa)
capture_90_ci = capture_radius_ci(disp_all, 0.9, radius,taxa= taxa)

ci_tot <- cbind(capture_50_ci,capture_70_ci[,3:5],capture_90_ci[,3:5])

dispersal_distance[[ii-2]] <- ci_tot
}

write.csv(dispersal_distance[[1]],paste(data.loc,'catchment_radii_var_Gaussian.csv',sep=''))
write.csv(dispersal_distance[[2]],paste(data.loc,'catchment_radii_var_PLK.csv',sep=''))

saveRDS(dispersal_distance,paste(data.loc,'catchment_radii.RDS',sep=''))















# 
# alnus.index <- grep('Alnus',colnames(pollen.test))
# colSums(pollen.test[,alnus.index])
# alnus.tot <- rowSums(pollen.test[,alnus.index])
# 
# corylus.index <- grep('Corylus',colnames(pollen.test))
# corylus.tot <- rowSums(pollen.test[,corylus.index])
# 
# ostrya.index <- grep('Ostrya',colnames(pollen.test))
# ostrya.tot <- rowSums(pollen.test[,ostrya.index])
# 
# salix.index <- grep('Salix',colnames(pollen.test))
# salix.tot <- (pollen.test[,salix.index])
# 
# tilia.index <- grep('Tilia',colnames(pollen.test))
# tilia.tot <- rowSums(pollen.test[,tilia.index])
# 
# 
# pinus.index <- grep('Pinus',colnames(pollen.test))
# colSums(pollen.test[,pinus.index])
# pinus.tot <- rowSums(pollen.test[,pinus.index])
# 
# 
# load('~/workflow_stepps_calibration/calibration/data/elicitation_neus_certainty_median_79_sites_only_abies.RData')
# cbind(calib_trans$`Other hardwood`[order(calib_trans$dataset.id)],alnus.tot)
# cbind(calib_trans$`Other hardwood`[order(calib_trans$dataset.id)],corylus.tot)
# cbind(calib_trans$`Other hardwood`[order(calib_trans$dataset.id)],ostrya.tot)
# 
# test <- cbind(calib_trans$`Other hardwood`[order(calib_trans$dataset.id)],rowSums(cbind(alnus.tot,ostrya.tot)))
# ratio <- round(test[,2]/test[,1],2)
# hist(ratio)
# sum(ratio>=0.6)
# 
# 
# cbind(sort(y[,"Pine"]),sort(pinus.tot))
# 
# cbind
# 
# oh_differentiated <- matrix(ncol=length(other_hardwood),nrow=nrow(pollen.test))
# for(i in other_hardwood){
#   index <- grep(i,colnames(pollen.test))
#   if(length(index)>1) {oh_differentiated[,other_hardwood==i] <- rowSums(pollen.test[,index])}
#   if(length(index)==1) {oh_differentiated[,other_hardwood==i] <- pollen.test[,index]}
# }
# 
# colnames(oh_differentiated) <- other_hardwood
# cs <-   colSums(oh_differentiated)
# 
# 
# cbind(calib_trans$`Other hardwood`[order(calib_trans$dataset.id)],rowSums(oh_differentiated))
