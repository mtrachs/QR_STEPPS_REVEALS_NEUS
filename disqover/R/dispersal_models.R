##################################################################################
##### DispersalFaktorK
##################################################################################

#' @export
# calculates the dispersal factor K for each species (=influx if the cover of that taxon were 100%)
DispersalFactorK=function(vg, tBasin, dBasin, dwm, regionCutoff, u = 3){

  rBasin <- dBasin/2

  # select chosen dispersal model (for parameters of GPM model see Jackson & Lyford 1999)
  if (dwm=="gpm neutral") {
    disModel<-gpm(fallSpeed=vg, cz=0.12, n=0.25, u = u)
  }
  else if (dwm=="gpm unstable") {
    disModel<-gpm(fallSpeed=vg, cz=0.21, n=0.20, u = u)
  }
  else if (dwm=="1overd") {
    disModel<-OneOverD.f()
  }
  else if (dwm=="lsm unstable")   {
    vgr<-100*round(vg,2)+1 #find right column
    disModel<-lsm_unstable[[vgr]]

  # with the following option, alternative models may be used...
  #} else if (dwm=="alternative")   {
  #  load("alternative.rda") # load model, which is prepared using a loess smoother
  #  vgr<-100*round(vg,2)+1 #find right column
  #  disModel<-alternative[[vgr]]

  }
  else stop("no valid distance weighting method selected")


  # for lakes, add additinal deposition arriving due to lake mixing
  if (tBasin == "peatland") {

    # deposition in a peatland
    if (regionCutoff <= 1e6) influx <- predict(disModel,rBasin)-predict(disModel,regionCutoff)

    # if regionCutOff > 1000 km - ignore pollen from beyond that limit
    else influx <- predict(disModel,rBasin)
  }
  else if (tBasin == "lake") {

    # set ring width for lake model, 2 m in smaller basins and 10 m in larger basins
    stepSize=2
    if (dBasin>=2000) stepSize=10

    # create vector with steps along lake radius from center to margin
    rSteps <- seq(from = stepSize, to = rBasin-1, by = stepSize)

    # call the lake model
    lakeInflux <- do.call(rbind,lapply(rSteps, LakeModel, dBasin=dBasin, disModel=disModel, stepSize=stepSize, regionCutoff=regionCutoff))

    #devide total influx by (relevant) lake area
    influx <- sum(lakeInflux) / (pi*(max(rSteps))^2)
  }

  return(influx)
}


##################################################
##### LakeModel
##################################################

# calculates additonal pollen deposition on lake points between the center and margin of the lake
# x is the distance between lake center and point for which deposition is calculated

LakeModel = function(x, dBasin, disModel, stepSize, regionCutoff){
  rBasin<-dBasin/2

  # create a sequence of radii from the smallest to the largest needed
  r2 <- seq(from = rBasin - x, to = 2*rBasin, by = 2)

  # cut off high values (>rBasin+x), otherwise alpha calculation fails
  r2Cut <- r2; r2Cut[r2Cut>(rBasin+x)]=rBasin+x

  # calculate the angel alpha (see Appendix)
  alpha <- acos((r2Cut^2+x^2-rBasin^2)/(2*r2Cut*x))

  # the proportion of the circle outside the lake area
  propOut <- 1-alpha/pi
  a <- diff(propOut)/2 ; propOut <- propOut[-1]-a   # produce mean values between two circles

  # pollen airborne at each circle
  airborne <- predict(disModel, r2)

  # calculate difference -> influx from full rings
  influxRing <- abs(diff(airborne))

  # multiply with proportion of each ring that is outside the lake --> true influx from each ring
  influxRing <- influxRing*propOut

  # influx is sum of all rings plus deposition from outside
  influx <- sum(influxRing) + predict(disModel,2*rBasin) - predict(disModel,regionCutoff)

  ringArea <- pi*x^2-pi*(x-stepSize)^2

  # return total influx multiplied by ringArea to account for the increasing area with increasing x
  return(influx*ringArea)
}


##########################################################################
##### Gaussian plume model (GPM) - Sutton model follwing Prentice (1985)
##########################################################################

# prepares distance weighting factors with a Gaussian plume model
# Arguments:
# n: turbulence parameter
# cz: vertical diffusion coefficient
# u: mean wind speed [meter/s]

gpm = function (fallSpeed, cz, n, u){
  x <- seq(from = 0, to = 100, by = 0.1)
  x <- x^3
  y <- exp((-4*fallSpeed*x^(n/2))/(n*u*sqrt(pi)*cz))
  return(loess(y ~ x, span=0.015, model=FALSE, degree=2))
}


###################################################
##### 1overd dispersal model
###################################################

# prepares a model with the 1/d dispersal function
# bitte ignorieren - noch nicht fertig
OneOverD.f = function (){
  x <- seq(from = 1, to = 100, by = 0.15) # klappt nicht sauber
  x <- x^3
  y <- 1/x
  # es muss ein remaining airborne vector bei raus kommen
  y<-y/sum(y)
  # plot(y~x, xlim=c(0,100))
  loess(y ~ x, span=0.015)
  return(loess(y ~ x, span=0.015, model=FALSE, degree=2))
}
