#--------------------------------------------------------------------------------------------------------------------
library(rstan)
library(fields)
library(abind)

setwd(paste0(wd,'figures/'))
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'figures/'





pl_kernel <- function(d,a,b) {(b-2)*(b-1)/(2*pi*a^2)*(1+d/a)^-b}
g_kernel <- function(d,psi) {exp(-(d^2/psi^2))}

taxa <- c('Ash','Beech','Birch','Elm','Hemlock','Maple','Oak','Other conifer',
               'Other hardwood','Pine','Spruce','Tamarack')


#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#generate d_pot

locations <- seq(8,700,8)
loc.tot <- c(rev(-locations),0,locations)
stepps.grid <- expand.grid(loc.tot,loc.tot)

distances <- rdist(matrix(c(0,0),ncol=2),stepps.grid)
dist.num <- table(distances)

column1 <- as.numeric(names(dist.num))

d_pot <- cbind(column1,unname(dist.num))
d_pot <- d_pot[((d_pot[,1]>0)&(d_pot[,1]<700)),]
d_pot[,1] <- d_pot[,1]/1000
d_pot <- rbind(d_pot,c(0.8,0))
#--------------------------------------------------------------------------------------------------------------------
# Figure comparing  UMW and NEUS
#--------------------------------------------------------------------------------------------------------------------
taxa_umw <- c('Ash','Beech','Birch','Elm','Hemlock','Maple','Oak','Other conifer',
              'Other hardwood','Pine','Spruce','Tamarack')


taxa_neus <- c('Ash','Beech','Birch','Chestnut','Hemlock','Hickory','Maple','Oak','Other conifer',
               'Other hardwood','Pine','Spruce','Tamarack')

taxa_shared <- intersect(taxa_umw,taxa_neus)

source('R/build_cal_main.r') # 

taxa.important <- c('Beech','Hemlock','Oak','Spruce')

pdf(paste(plot.loc,'Figure_7.pdf',sep=''),height=8,width=8)
#pdf(paste(plot.loc,'Figure_7_all_taxa.pdf',sep=''),height=16,width=12)
#par(mfrow=c(4,3))
par(mfrow=c(2,2),mar=c(1,0,1,0),oma=c(3,4,1,1))
#for(taxon in taxa_shared) {  
for(taxon in taxa.important) {  
  
  
  plot(c(1,1),type='n',xlim=c(0,800),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  mtext(side=3,line=0,font=2.2,text=taxon)
  box()
  if(taxon==taxa.important[1]){
      axis(2)
    
    mtext(side=2,line=2.2, text = 'Proportion of deposited pollen',font = 2,cex = 1)
    }
  if(taxon==taxa.important[3]){
    axis(1)
    axis(2)
    mtext(side=1,line=2.2, text = 'distance [km]',font = 2,cex = 1)
    mtext(side=2,line=2.2, text = 'Proportion of deposited pollen',font = 2,cex = 1)
  }
  if(taxon==taxa.important[4]){
    axis(1)
    mtext(side=1,line=2.2, text = 'distance [km]',font = 2,cex = 1)
  }
  if(taxon=='Beech'){
  #legend('bottomright',col=c(1,1,rep('green',2),rep(4,2),rep(2,2)),lty=rep(c(1,2),4),
  #       legend = paste(rep(c('Base GK','Variable GK','Base PLK','Variable PLK'),2),
 #                     rep(c('UMW','NEUS'),each=4)),cex = 0.95,lwd = 2)
    legend('bottomright',col=rep(c(1,2),each=2),lty=rep(c(1,2),2),
           legend = paste(rep(c('Variable PLK','Variable GK'),2),
                         rep(c('UMW','NEUS'),each=2)),cex = 0.95,lwd = 2)
  }
    
    
    
  for(region in c('UMW','NEUS')){

runs1 <- list(runs[[3]],runs[[4]])

for (run in runs1) {
  kernel <- run$kernel
  num_a <- run$one_a
  one_psi <- run$one_psi
  handle <- strsplit(run$suff_fit,'_A')[[1]][1]
  suff_fit <- run$suff_fit
  
  if(region=='UMW'){
  fname <- read_stan_csv(paste(data.loc,suff_fit,'.csv',sep=''))
  post <- rstan::extract(fname,permuted=FALSE)
  }
  
  if(region=='NEUS') {
    
  #look at that again...
    # would have to change this if not run_pl
      if(kernel=='pl'){
          fname <- paste0(wd,'results/data/cal_pl_Ka_Kgamma_EPs_modified_a_110.csv')
      }
      if(kernel=='gaussian'){
          fname <- paste0(wd,'results/data/cal_g_Kpsi_Kgamma_EPs_110.csv')
        }
      fit <- read_stan_csv(fname)
      post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
      }
    

  
  
  #look at that again...
  param.names <-colnames(post[,1,]) #find parameter names of posterior
  param.kernel <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
  gamma <- as.matrix(post[,1,param.kernel%in%'gamma'])
  
  if((ncol(gamma)>1)&(region=='NEUS')){
    gamma <- gamma[,taxa_neus%in%taxa_umw]
  }
  
  if(ncol(gamma)>1){
    gamma <-colMeans(gamma)
  } else{
  gamma <- mean(gamma)
  }
  
  
  
  if(kernel =='pl') {
    param.kernel1 <- param.kernel%in%c('a','b') 
    est.kernel <- apply(post[,1,param.kernel1],2,median)
    a <- est.kernel[grep('a',names(est.kernel))]
    b <- est.kernel[grep('b',names(est.kernel))]
    
    if((region=='NEUS')&(length(a)>1)){a <- a[taxa_neus%in%taxa_umw]}
    
    if((region=='NEUS') & (length(a)>1)){y <- a[taxon==taxa_shared]}  
    if((region=='UMW') & (length(a)>1)){y <- a[taxon==taxa_umw]}  
    if(length(a)==1){y <- a}
  
      kernel_shape <- sapply(d_pot[,1],function(x){pl_kernel(x,y,b)})
      pollen.weight <- d_pot[,2] * pl_kernel(d_pot[,1],y,b)
      
      
      
      
      
      gamma1 <- gamma[a==y]
      
      points(1000*d_pot[,1],gamma1 + (1-gamma1)*cumsum(pollen.weight)/sum(pollen.weight),type='l',lwd=2,
          xlab='',ylab='',lty = 1,col=ifelse(region=='UMW',1,2))
  }
  
  if(kernel == 'gaussian') {
    param.kernel1 <- param.kernel%in%'psi'
    if(sum(param.kernel1)>1) {psi <- apply(post[,1,param.kernel1],2,median)}
    if(sum(param.kernel1)==1) {psi <- median(post[,1,param.kernel1])}
    
    if((region =='NEUS')&(length(psi)>1)){psi <- psi[taxa_neus%in%taxa_umw]}
    
    if((region =='NEUS')& (length(psi)>1)){y <- psi[taxon == taxa_shared]}  
    if((region=='UMW') & (length(psi)>1)){y <- psi[taxon == taxa_umw]} 
    if(length(psi)==1){y <- psi}

      kernel_shape <- sapply(d_pot[,1],function(x){g_kernel(x,y)})
      pollen.weight <- d_pot[,2] * g_kernel(d_pot[,1],y)
      
      gamma1 <- gamma[psi==y]
      
      points(1000*d_pot[,1],gamma1 + (1-gamma1)*cumsum(pollen.weight)/sum(pollen.weight),type='l',lwd=2,
          xlab='',ylab='',lty = 2,col=ifelse(region=='UMW',1,2))

  }
  }
  }
}
dev.off()



##########################################################################################################################
# Same for all taxa
##########################################################################################################################
pdf(paste(plot.loc,'SF6.pdf',sep=''),height=16,width=12)
par(mfrow=c(4,3),mar=c(1,0,1,0),oma=c(3,4.25,1,1),cex.axis= 2)
for(taxon in taxa_shared) {  

  
  plot(c(1,1),type='n',xlim=c(0,800),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
  mtext(side=3,line=0,font=2.2,text=taxon, cex = 1.5)
  box()
  if(taxon%in%taxa_shared[c(1,4,7)]){
    axis(2)
    mtext(side=2,line=2.6, text = 'Proportion of deposited pollen',font = 2,cex = 1.3)
  }
  if(taxon%in%taxa_shared[10]){
    axis(1)
    axis(2)
    mtext(side=1,line=2.5, text = 'distance [km]',font = 2,cex = 1.5)
    mtext(side=2,line=2.6, text = 'Proportion of deposited pollen',font = 2,cex = 1.3)
  }
  if(taxon%in%taxa_shared[c(9,10,11)]){
    axis(1)
    mtext(side=1,line=2.5, text = 'distance [km]',font = 2,cex = 1.5)
  }
  if(taxon=='Beech'){
    #legend('bottomright',col=c(1,1,rep('green',2),rep(4,2),rep(2,2)),lty=rep(c(1,2),4),
    #       legend = paste(rep(c('Base GK','Variable GK','Base PLK','Variable PLK'),2),
    #                     rep(c('UMW','NEUS'),each=4)),cex = 0.95,lwd = 2)
    legend('bottomright',col=rep(c(1,2),each=2),lty=rep(c(1,2),2),
           legend = paste(rep(c('Variable PLK','Variable GK'),2),
                          rep(c('UMW','NEUS'),each=2)),cex = 0.95,lwd = 2)
  }
  
  
  
  for(region in c('UMW','NEUS')){
    
    runs1 <- list(runs[[3]],runs[[4]])
    
    for (run in runs1) {
      kernel <- run$kernel
      num_a <- run$one_a
      one_psi <- run$one_psi
      handle <- strsplit(run$suff_fit,'_A')[[1]][1]
      suff_fit <- run$suff_fit
      
      if(region=='UMW'){
        fname <- read_stan_csv(paste(data.loc,suff_fit,'.csv',sep=''))
        post <- rstan::extract(fname,permuted=FALSE)
      }
      
      if(region=='NEUS') {
        
        #look at that again...
        # would have to change this if not run_pl
        if(kernel=='pl'){
          fname <- paste0(wd,'results/data/cal_pl_Ka_Kgamma_EPs_modified_a_110.csv')
        }
        if(kernel=='gaussian'){
          fname <- paste0(wd,'results/data/cal_g_Kpsi_Kgamma_EPs_110.csv')
        }
        fit <- read_stan_csv(fname)
        post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
      }
      
      
      
      
      #look at that again...
      param.names <-colnames(post[,1,]) #find parameter names of posterior
      param.kernel <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
      gamma <- as.matrix(post[,1,param.kernel%in%'gamma'])
      
      if((ncol(gamma)>1)&(region=='NEUS')){
        gamma <- gamma[,taxa_neus%in%taxa_umw]
      }
      
      if(ncol(gamma)>1){
        gamma <-colMeans(gamma)
      } else{
        gamma <- mean(gamma)
      }
      
      
      
      if(kernel =='pl') {
        param.kernel1 <- param.kernel%in%c('a','b') 
        est.kernel <- apply(post[,1,param.kernel1],2,median)
        a <- est.kernel[grep('a',names(est.kernel))]
        b <- est.kernel[grep('b',names(est.kernel))]
        
        if((region=='NEUS')&(length(a)>1)){a <- a[taxa_neus%in%taxa_umw]}
        
        if((region=='NEUS') & (length(a)>1)){y <- a[taxon==taxa_shared]}  
        if((region=='UMW') & (length(a)>1)){y <- a[taxon==taxa_umw]}  
        if(length(a)==1){y <- a}
        
        kernel_shape <- sapply(d_pot[,1],function(x){pl_kernel(x,y,b)})
        pollen.weight <- d_pot[,2] * pl_kernel(d_pot[,1],y,b)
        
        
        
        
        
        gamma1 <- gamma[a==y]
        
        points(1000*d_pot[,1],gamma1 + (1-gamma1)*cumsum(pollen.weight)/sum(pollen.weight),type='l',lwd=2,
               xlab='',ylab='',lty = 1,col=ifelse(region=='UMW',1,2))
      }
      
      if(kernel == 'gaussian') {
        param.kernel1 <- param.kernel%in%'psi'
        if(sum(param.kernel1)>1) {psi <- apply(post[,1,param.kernel1],2,median)}
        if(sum(param.kernel1)==1) {psi <- median(post[,1,param.kernel1])}
        
        if((region =='NEUS')&(length(psi)>1)){psi <- psi[taxa_neus%in%taxa_umw]}
        
        if((region =='NEUS')& (length(psi)>1)){y <- psi[taxon == taxa_shared]}  
        if((region=='UMW') & (length(psi)>1)){y <- psi[taxon == taxa_umw]} 
        if(length(psi)==1){y <- psi}
        
        kernel_shape <- sapply(d_pot[,1],function(x){g_kernel(x,y)})
        pollen.weight <- d_pot[,2] * g_kernel(d_pot[,1],y)
        
        gamma1 <- gamma[psi==y]
        
        points(1000*d_pot[,1],gamma1 + (1-gamma1)*cumsum(pollen.weight)/sum(pollen.weight),type='l',lwd=2,
               xlab='',ylab='',lty = 2,col=ifelse(region=='UMW',1,2))
        
      }
    }
  }
}
dev.off()

#---------------------------------------------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------------------------------------------

# if(region=='NEUS')
# 
#   kernel <- run$kernel
#   num_a <- run$one_a
#   one_psi <- run$one_psi
#   handle <- strsplit(run$suff_fit,'_A')[[1]][1]
#   
#   #look at that again...
#   for (i in 0:9) { # would have to change this if not run_pl
#     if(kernel=='pl'){
#       if(num_a==FALSE){
#         fname = paste('~/stepps_median/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
#       }
#       if((num_a==TRUE)&(i<5)){
#         fname = paste('~/workflow_stepps_calibration/results/data/stepps_median/output_cal_pl_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
#       }
#       if(run$a_const==TRUE){
#         fname = paste('~/stepps_modified_stan/output/output_cal_pl_Ka_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
#       }  
#     }
#     if(kernel=='gaussian'){
#       if(one_psi==FALSE){
#         fname = paste('~/stepps_median/output_cal_g_Kpsi_Kgamma_EPs_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
#       }
#       if(one_psi==TRUE){
#         fname = paste('~/stepps_median/output_neus_',i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
#       }}
#     fit <- read_stan_csv(fname)
#     if (i==0)  {post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)}
#     else{post_new <-  rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
#     post <- abind(post,post_new,along=1)
#     }
#   }
# }
#   
#   if(kernel =='pl') {
#     param.names <-colnames(post[,1,]) #find parameter names of posterior
#     param.kernel <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
#     param.kernel1 <- param.kernel%in%c('a','b') 
#     est.kernel <- apply(post[,1,param.kernel1],2,median)
#     a <- est.kernel[grep('a',names(est.kernel))]
#     b <- est.kernel[grep('b',names(est.kernel))]
#     
#     sapply(a,function(y){
#       kernel_shape <- sapply(distances,function(x){pl_kernel(x,y,b)})
#       pollen.weight <- d_pot[,2] * pl_kernel(d_pot[,1],y,b)
#       plot(1000*distances,kernel_shape,type='b',pch = 16,cex = 0.75,main = paste (taxa[y==a],handle),
#            xlab='',ylab='')
#       mtext(side=1,line=2.2, text = 'distance',font = 2)
#       mtext(side=2,line=2.2, text = 'Proportion of deposited pollen',font = 2)
#       plot(1000*d_pot[,1],cumsum(pollen.weight)/sum(pollen.weight),type='l',lwd=2,
#            main = paste (taxa[y==a],handle),xlab='',ylab='')
#       mtext(side=1,line=2.2, text = 'distance',font = 2)
#       mtext(side=2,line=2.2, text = 'Proportion of deposited pollen',font = 2)
#     })
# 
#   }
#   
#   if(kernel =='gaussian') {
#     param.names <-colnames(post[,1,]) #find parameter names of posterior
#     param.kernel <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
#     param.kernel1 <- param.kernel%in%'psi'
#     if(sum(param.kernel1)>1) {psi <- apply(post[,1,param.kernel1],2,median)}
#     if(sum(param.kernel1)==1) {psi <- median(post[,1,param.kernel1])}
#     
#     sapply(psi,function(y){
#       kernel_shape <- sapply(distances,function(x){g_kernel(x,y)})
#       pollen.weight <- d_pot[,2] * g_kernel(d_pot[,1],y)
#       plot(1000*distances,kernel_shape,type='b',pch = 16,cex = 0.75,main = paste (taxa[y==psi],handle),
#            ylim=c(0,1),xlab='',ylab='')
#       mtext(side=1,line=2.2, text = 'distance',font = 2)
#       mtext(side=2,line=2.2, text = 'weight',font = 2)
#       plot(1000*d_pot[,1],cumsum(pollen.weight)/sum(pollen.weight),type='l',lwd=2,
#            main = paste (taxa[y==psi],handle),xlab='',ylab='')
#       mtext(side=1,line=2.2, text = 'distance',font = 2)
#       mtext(side=2,line=2.2, text = 'Proportion of deposited pollen',font = 2)
#     })
#   }
# }
# dev.off()  
