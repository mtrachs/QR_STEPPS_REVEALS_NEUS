library(maps)
plot.loc <- paste0(wd,'figures/figures/')


corners =matrix(ncol=2,byrow=TRUE, 
                c(-81,47.5,
                  -66.5,47.5,
                  -66.5,39.5,
                  -81,39.5,
                  -81,47.5))

outer.ticks1 <- matrix(ncol=2,byrow=FALSE,
                      c(rep(-170,6),seq(25,75,10)))

outer.ticks2 <- matrix(ncol=2,byrow=FALSE,
                       c(rep(-50,6),rev(seq(25,75,10))))

outer.ticks3 <- matrix(ncol=2,byrow=FALSE,
                       c(rev(c(seq(-50,-170,-25),-170)),rep(75,6)))

outer.ticks4 <- matrix(ncol=2,byrow=FALSE,
                       c(c(seq(-50,-170,-25),-170),rep(25,6)))

outer.ticks <- rbind(outer.ticks1,outer.ticks3,outer.ticks2,outer.ticks4)

  
pdf(paste0(plot.loc,'na_inset.pdf'),height=10,width = 16)
  par(mar=c(0,0,0,0))
  map(regions=c('canada','usa','mexico'),xlim=c(-180,-40),ylim=c(22,76),col='white')
  map(regions=c('canada','usa','mexico'),xlim=c(-180,-40),ylim=c(25.75,74),add=TRUE)
  points(corners,pch = 20,type='l',lwd=10,col='blue')
  text(-100,55,font=2,'Canada',pos=1,cex=3.5)
  text(-100,40,font=2,'United States',pos=1,cex=3.5)
  text(-108,29,font=2,'Mexico',pos=1,cex=3.5)
  points(outer.ticks,type='b',pch=3)
  text(outer.ticks2[1:6,],paste(as.character(abs(outer.ticks2[1:6,2])), "\u00b0","N", sep=''),pos=4,font=2,cex=3)
  text(outer.ticks4[1:6,],paste(as.character(abs(outer.ticks4[1:6,1])), "\u00b0","W", sep=''),pos=1,font=2,cex=3,offset=1)
dev.off()  
