trace_plots <- function(fit, pars, N_cores, suff="", save_plots=TRUE, fpath,fname){
  
  post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  
  n      = dim(post)[3]
  #idx = c(seq(1,n-N_cores-1), n)
  idx = pars
  labels = colnames(post[,1,idx])
  
  ntrace = length(idx)
  
  avg = summary(fit)$summary[,"50%"]
  
  if (save_plots){
    if (nchar(suff)>1){
      suff = paste0('_', suff)
    }
    pdf(paste(fpath, "/",fname, suff, ".pdf", sep=""), width=8, height=6)
  }
  
  par(mfrow=c(1,1))
  for (i in 1:ntrace){
    
    plot(post[,1,idx[i]], type="l", ylab=labels[i], xlab="iter")
    abline(h=avg[idx[i]], col="blue")
    abline(h=summary(fit)$summary[,"2.5%"][idx[i]], col='blue', lty=2)
    abline(h=summary(fit)$summary[,"97.5%"][idx[i]], col='blue', lty=2)
  }
  
  if (save_plots){
    dev.off()
  }
}
