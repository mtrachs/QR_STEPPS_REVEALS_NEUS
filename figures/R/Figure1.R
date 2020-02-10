# library(sp)
# library(rgdal)
library(maptools)

#-------------------------------------------------------------------------------------------------------------------
setwd(paste0(wd,'calibration/'))
help.fun.loc <- 'calibration_helper_funs/'
data.loc <- 'data/'
plot.loc <- paste0(wd,'figures/figures/')
source('utils/dataPlotFuns.r')
source('utils/simDataFuns.r')

#have to change that one
load(paste0(data.loc,'/elicitation_neus_certainty_median_110_sites_only_abies_new_species.RData'))

########################################################################################################
## pls pie map
########################################################################################################
# reorder data by pollen abundance
new.order = order(colSums(y), decreasing=TRUE)
taxa <- colnames(y)

taxa = taxa[new.order]
y = y[,new.order]


r = r[,new.order]
colnames(r) = taxa

col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
             "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",'cyan',"#CAB2D6")

stepps_cols = data.frame(taxa=taxa, cols=col_list)

########################################################################################################
## pls pie map
########################################################################################################

ntaxa = length(taxa)
#what is this
centers   = veg_coords#centers_veg
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])

subgrid = regular_subgrid(centers, dx=34000,dy=34000)
knots_in = knots_in_domain4(subgrid, centers, cell_width = 8000)

plot(centers[,1], centers[,2], asp=1)
points(knots_in[,1], knots_in[,2], col='red')

d = rdist(knots_in, centers)

pls_coarse = matrix(0, nrow=nrow(d), ncol=ntaxa)
colnames(pls_coarse) = taxa

for (i in 1: nrow(centers)){
  #print(min(d[,i]))
  close_knot = which.min(d[,i])
  
  #   print(as.double(pls_cut[i,]))
  #   print(pls_cut[i,])
  
  
  pls_coarse[close_knot,] <- pls_coarse[close_knot,] + as.matrix(r[i,])
}

pls_props  = t(apply(pls_coarse, 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))

shift=30000


#setwd('C:/Work/Paleon/STEPPS/stepps-calibration-master/')
# postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
pdf(paste(plot.loc,'Fig1b.pdf',sep=''), width=12, height=10)
par(mfrow=c(1,1))
pieMap(proportions = pls_props, 
       centers  = knots_in,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 14000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y',
       add_legend = TRUE, 
       main_title='',
       col_list= col_list,
       cex.legend = 1.7,
       state.plot.index = 1:20)
dev.off()


########################################################################################################
## pollen pie map
########################################################################################################

pollen = y
pollen_props  = t(apply(pollen, 1, function(x) x/sum(x)))
colnames(pollen_props) = taxa


sputm <- SpatialPoints(pollen_coords, proj4string=CRS('+init=epsg:4326'))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))


centers   = spgeo@coords#centers_polA
colnames(centers) = c('x', 'y')


#############################################################################################################
# coordinates of bounding boxes 
#############################################################################################################
yb <- seq(39.5,47.5,0.1)
xb <- seq(-81,-66.5,0.1)


w_bound <- cbind(rep(-81,length(yb)),yb) 
n_bound <- cbind(xb,rep(47.5,length(xb))) 
e_bound <- cbind(rep(-66.5,length(yb)),rev(yb)) 
s_bound <- cbind(rev(xb),rep(39.5,length(xb))) 

bound <- rbind(w_bound,
               n_bound,
               e_bound,
               s_bound)

sputm <- SpatialPoints(bound, proj4string=CRS('+init=epsg:4326'))
bounds_us <- spTransform(sputm, CRS("+init=epsg:3175"))


bounding_corners <- matrix(ncol=2,c(-81,39.5,
                                    -81,47.5,
                                    -66.5,47.5,
                                    -66.5,39.5),
                           byrow = TRUE)

sputm <- SpatialPoints(bounding_corners, proj4string=CRS('+init=epsg:4326'))
bounding_corners_us <- spTransform(sputm, CRS("+init=epsg:3175"))


bounds_us <- bounds_us@coords 
bounding_corners_us <- bounding_corners_us@coords

#####################################################################################################################

pdf(paste(plot.loc,'Fig1a.pdf',sep=''), width=12, height=10)
par(mfrow=c(1,1))
pieMap(proportions = pollen_props, 
       centers  = centers,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo-shift, xhi+3*shift),
       ylim   = c(ylo-shift, yhi+2.5*shift),
       radius = 18000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y', 
       add_legend=FALSE,
       main_title='',
       col_list=col_list,
       state.plot.index = 1:20)
points(bounds_us,pch = 20,type='l',lwd=1.5)
text(bounding_corners_us[1,1],bounding_corners_us[1,2],labels='81W',pos=1,cex=1)
text(bounding_corners_us[1,1],bounding_corners_us[1,2],labels='39.5N',pos=2,cex=1)
text(bounding_corners_us[2,1],bounding_corners_us[2,2],labels='81W',pos=3,cex=1)
text(bounding_corners_us[2,1],bounding_corners_us[2,2],labels='47.5N',pos=2,cex=1)
text(bounding_corners_us[3,1],bounding_corners_us[3,2],labels='66.5W',pos=3,cex=1)
text(bounding_corners_us[3,1],bounding_corners_us[3,2],labels='47.5N',pos=4,cex=1)
text(bounding_corners_us[4,1],bounding_corners_us[4,2],labels='66.5W',pos=1,cex=1)
text(bounding_corners_us[4,1],bounding_corners_us[4,2],labels='39.5N',pos=4,cex=1)
dev.off()
# pdf('C:/Work/Paleon/STEPPS/NEUS/pie_plot_pollen_NEUS_ensemble_new.pdf', width=12, height=10)
# 
#   sapply(1:100,function(zzz) {
# 
#     load(paste('C:/Work/Paleon/STEPPS/NEUS/data_for_stan/elicitation_neus_',zzz,'.RData',sep=''))
# 
#     y = y[,new.order]
#     pollen = y
#     pollen_props  = t(apply(pollen, 1, function(x) x/sum(x)))
#     colnames(pollen_props) = taxa
#     sputm <- SpatialPoints(pollen_coords, proj4string=CRS('+init=epsg:4326'))
#     spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
#     centers   = spgeo@coords#centers_polA
#     colnames(centers) = c('x', 'y')
# 
# 
# 
#   par(mfrow=c(1,1))
#   pieMap(proportions = pollen_props, 
#        centers  = centers,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 18000,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y', 
#        add_legend=TRUE,
#        main_title='',
#        col_list=col_list)
#   })
# dev.off()
