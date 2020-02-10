dispersal_cdf <- function(radius, run, taxa, path_out,d_pot,post,rescale){
  
  K = length(taxa)  
  
  #coord_pot = seq(-700000, 700000, by=8000)
  #coord_pot = expand.grid(coord_pot, coord_pot)
  
  #d_pot = t(rdist(matrix(c(0,0), ncol=2), as.matrix(coord_pot, ncol=2))/rescale)
  
  # want a circular region
  #idx_circ  = which(d_pot[,1] < 0.7)
  #coord_pot = coord_pot[idx_circ, ]
  #d_pot     = d_pot[idx_circ, ]
  
  #d_pot = unname(as.matrix(count(d_pot)))
  N_pot = nrow(d_pot)
  
  dat = data.frame(matrix(0, nrow=0, ncol=4))
  r_all = list(length=length(runs))
  #perhaps have run and post outside 
    
    #     sum_w = build_sumw_pFot(post, K, N_pot, d_pot, run)
    #     r_int = dispersal_decay(post, d_pot, sum_w, radius/rescale, run, taxa)
    #     colnames(r_int) = taxa
    
    sum_w = build_sumw_pot_ci(post, K, N_pot, d_pot, run)
    r_niter = dispersal_decay_ci(post, d_pot, sum_w, radius/rescale, run, taxa)
    
    r_all = list(r_niter, run$handle)
    #     r_all[[i]]$handle = run$handle
    
    r_int = t(apply(r_niter, 1, rowMeans))
    colnames(r_int) = taxa
    
    run_dat = data.frame(radius = radius/1e3, r_int, handle=rep(run$handle, length(radius)))
    dat = rbind(dat, melt(run_dat, id=c('radius', 'handle')))
  
  
  levels(dat$handle) <- "Variable PL" #c("Base G", "Base PL", "Variable G", "Variable PL")
  
  return(list(dat=dat, r_all=r_all))
}

plot_dispersal_cdfs <- function(dat,fig_name,x.lim){
  
  levels(dat$variable)[levels(dat$variable) == "OTHER.CONIFER"] = "FIR"
  levels(dat$variable)[levels(dat$variable) == "OTHER.HARDWOOD"] = "OTHER HARDWOOD"
  
  dat$variable = factor(dat$variable, levels = sort(levels(dat$variable)))
  
  p <- ggplot(dat) + geom_line(data=dat, 
                               aes(x=radius, y=value, colour=factor(handle), linetype=factor(handle)), 
                               lwd=1.0) 
  p <- p + scale_colour_manual("Model", values=c('gray5', 'grey57', 'gray5', 'grey57')) 
  p <- p + scale_linetype_manual("Model", values = c("solid", "solid", "dashed", "dashed")) +
    xlab('Radius (km)') + ylab('Proportion of deposited pollen')
  p <- p + xlim(0, x.lim) + ylim(0,1.0)
  p <- p + facet_wrap(~variable, ncol=3) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                              panel.grid.minor = element_blank())
  print(p)
  
  ggsave(p, filename=sprintf('%s/%s/%s.pdf', wd, path_figs, fig_name), width=9, height=12)
  
}


capture_radius <- function(dat, prop){
  
  handles = as.vector(unique(dat$handle))
  
  cr = data.frame(handle=character(0), taxon=character(0), r=numeric(0))
  for (handle in handles){
    for (taxon in unique(dat$variable)){
      dat_sub = dat[which((dat$handle == handle) & (dat$variable == taxon)), ]
      rad = dat_sub[which.min(abs(dat_sub[, 'value'] - prop)), 'radius']
      cr = rbind(cr, data.frame(handle=handle, taxon=taxon, r=rad))
    }
  }
  
  return(cr)
}


capture_radius_ci <- function(dat, prop, radius,taxa){
  
  #   handles   = as.vector(unique(dat$handle))
  taxa_list = taxa
  ntaxa = length(taxa)
  
  cr = data.frame(handle=character(0), taxon=character(0), r_med=numeric(0), r_lb=numeric(0), r_ub=numeric(0))
  for (i in 1:(length(dat)-1)){
    for (j in 1:ntaxa){
      dat_sub = dat[[1]][,j,]
      cr_all    = apply(dat_sub, 2, function(x) radius[which.min(abs(x - prop))])
      cr_quants = quantile(cr_all, probs=c(0.025, 0.5, 0.975)) * 1e-3
      cr = rbind(cr, data.frame(handle = 'Variable PL',#'dat[[i]][[2]] 
                                taxon  = taxa_list[j], 
                                r_med  = cr_quants[2], 
                                r_lb   = cr_quants[1], 
                                r_ub   = cr_quants[3]))
    }
  }
  
  return(cr)
}


# do it!
