#####################################################################################################
# Figure to compare STEPPS and REVEALS PPEs and catchment radii and falling speeds
#####################################################################################################
library(abind)
library(rstan)

repo_location <- wd
setwd(repo_location)
stepps.data.loc <- 'results/data/'
reveasl.data.loc <- ''

catchment.radii <-readRDS(paste(stepps.data.loc,'catchment_radii.RDS',sep=''))

stepps.radii.plk <- catchment.radii[[2]]

####################################################################################################
#read REVEALS estimates
####################################################################################################
svs = read.csv(paste0(repo_location,'Reveals_NEUS/data/svs_LC6K.csv'), sep=',', header=TRUE, stringsAsFactors=FALSE)

svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)
svs_agg$taxon[svs_agg$taxon=='Alder'] <- 'Other hardwood'
svs_agg$taxon[svs_agg$taxon=='Fir'] <- 'Other conifer'
svs_agg$taxon[svs_agg$taxon=='Larch'] <- 'Tamarack'
svs_agg <- svs_agg[svs_agg$taxon%in%stepps.radii.plk$taxon,]


ppes <- readRDS(paste0(repo_location,'Reveals_NEUS/data/PPEs_agg.RDS'))
ppes$taxon <- as.character(ppes$taxon)
ppes$taxon[ppes$taxon=='Alder'] <- 'Other hardwood'
ppes$taxon[ppes$taxon=='Fir'] <- 'Other conifer'
ppes$taxon[ppes$taxon=='Larch'] <- 'Tamarack'
ppes <- ppes[ppes$taxon%in%stepps.radii.plk$taxon,]

stepps.radii.plk <- stepps.radii.plk[stepps.radii.plk$taxon%in%ppes$taxon,]

######################################################################################################
# load STEPPS pollen productivity estimates
######################################################################################################
fname <- paste0(repo_location,'results/data/cal_pl_Ka_Kgamma_EPs_modified_a_110.csv')
fit <- read_stan_csv(fname)
post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)

param.names <-colnames(post[,1,]) #find parameter names of posterior
param.greek <- sapply(param.names,function(x) strsplit(x,'[[]')[[1]][1]) #only take greek letter
param.e.bayes <- param.greek%in%c('gamma','phi')#only take eta and rho
remove <- grep('4',param.names[param.e.bayes]) # this is Chestnut
param.e.bayes[remove] <- FALSE
est.e.bayes <- colMeans(post[,1,param.e.bayes])
upper.e.bayes <- apply(post[,1,param.e.bayes],2,function(x) quantile(x,probs=0.975))
lower.e.bayes <- apply(post[,1,param.e.bayes],2,function(x) quantile(x,probs=0.0255))
phi <- est.e.bayes[grep('phi',names(est.e.bayes))]
upper.phi <- upper.e.bayes[grep('phi',names(upper.e.bayes))]
lower.phi <- lower.e.bayes[grep('phi',names(lower.e.bayes))]

names(phi) <- names(upper.phi) <- names(lower.phi) <- stepps.radii.plk$taxon 
phi.rescale <- phi/phi[names(phi)=='Oak'] 
upper.phi.rescale <- upper.phi/phi[names(phi)=='Oak'] 
lower.phi.rescale <- lower.phi/phi[names(phi)=='Oak'] 


#estimate ratio of Phis
phi.ensemble <- post[,1,param.greek=='phi'] 
phi.ratio <- phi.ensemble/phi.ensemble[,8]

names(phi) <- names(upper.phi) <- names(lower.phi) <- stepps.radii.plk$taxon 
phi.rescale <- colMeans(phi.ratio)[-4] #remove Chestnut
upper.phi.rescale <- apply(phi.ratio,2,function(x) quantile(x,probs=0.975))[-4] #remove Chestnut
lower.phi.rescale <- apply(phi.ratio,2,function(x) quantile(x,probs=0.025))[-4] #remove Chestnut

######################################################################################################
#load values of pollen dispersal depopsition coefficient (K)
K <- readRDS(paste0(wd,'Reveals_NEUS/data/K.RDS')) # have to change this
K <- round(K,4)
K <- as.data.frame(K)



######################################################################################################
# make plots two plots
######################################################################################################
setwd(paste(repo_location,'/figures',sep=''))
plot.loc <- 'figures/'



pdf(paste(plot.loc,'Figure5.pdf',sep=''),height = 10,width = 11)
par(mfrow=c(2,2))
par(oma=c(0,0,0,0.5),mar=c(3.5,3.5,1,0),cex.axis = 1.1) 

max.ppe <- max(c(phi.rescale,ppes$ppe[order(ppes$taxon)]))
plot(phi.rescale,ppes$ppe[order(ppes$taxon)],pch = 1:12,xlab='',ylab='',xlim=c(0,max.ppe ),ylim=c(0,max.ppe ),cex=1.5)
abline(a = 0,b= 1,lty = 2)
mtext(side=1,line=2.2,'Phi STEPPS',font=2)
mtext(side=2,line=2.2,'PPE REVEALS',font=2)
legend('topleft',legend = 'a)',box.col = NA,pt.cex=0,cex=1.5)



plot(stepps.radii.plk[,6],svs_agg$sv[order(svs_agg$taxon)],pch = 1:12,ylim=c(0,0.15),xlab='',ylab='',cex=1.5)
legend('topright',legend = stepps.radii.plk$taxon,pch = 1:12,cex=1.15)
mtext(side=1,line=2.2,'STEPPS 70% capture radius [km]',font=2)
mtext(side=2,line=2.2,'Fallspeed [m/s] ',font=2)
legend('topleft',legend = 'b)',box.col = NA,pt.cex=0,cex=1.5)

plot(stepps.radii.plk[,6],K$`550`,pch = 1:12,xlab='',ylab='',cex=1.5)
mtext(side=1,line=2.2,'STEPPS 70% capture radius [km]',font=2)
mtext(side=2,line=2.2,'Pollen dispersal-deposition coefficient',font=2)
mtext(side=3,line=0,'Lake radius: 550 m',font=2)
legend('topleft',legend = 'c)',box.col = NA,pt.cex=0,cex=1.5)

plot(K$`550`,K$`5.5`,pch = 1:12,xlab='',ylab='',ylim=c(0,max(K$`5.5`)),xlim=c(0,0.265),cex=1.5)
points(K$`550`,K$`55`,pch = 1:12,xlab='',ylab='',col=2,cex=1.5)
points(K$`550`,K$`5500`,pch = 1:12,xlab='',ylab='',col=4,cex=1.5)
mtext(side=1,line=2.2,'Pollen dispersal-deposition coefficient (550 m radius)',font=2)
mtext(side=2,line=2.2,'Pollen dispersal-deposition coefficient',font=2)
legend('topleft',legend = c(as.character(stepps.radii.plk$taxon),'5.5 m','55 m','5500'),
       ncol=2,pch = c(1:12,rep(16,3)),cex=1.2,col=c(rep(1,13),2,4))
legend('topright',legend = 'd)',box.col = NA,pt.cex=0,cex=1.5)


dev.off()





