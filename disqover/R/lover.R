
# am besten müssen die Daten wie bei MarcoPolo übergeben werden
# local und regional Pollen, assignment und parameter
#' @export


love <- function(loc, reg, params, dBasin, rsap, repeats = 100, dwm, tBasin = "peatland", nrsample){

  locsample <- unlist(loc[[nrsample+1]]) #a single local record
  regsample <- unlist(reg[[nrsample+1]]) #the corresponding regional record
  pro <- params[["PPEs"]] # pollen productivity/R-Values..



  #########################################
  # LOVE in Shinyas way
  #########################################

  love1 <- replicate(n=repeats,{

    # add random error to PPEs (using mean and stanadard deviation), 'while loop' is to avoid negative values
    PPEs <- rnorm_reveals(length(params[["PPEs"]]), params[["PPEs"]], params[["PPE.errors"]])

    #calculate the s roof from Sugita 2007 LOVE formula 7
    sr <- sr(regsample, params, tBasin = "peatland", dBasin=dBasin, rsap=rsap, dwm=dwm)

    # randomize local counts
    locsampler <- rmultinom(1, sum(locsample), locsample/sum(locsample))

    i <- locsampler/pro
    irel <- i/sum(i)

    fin <- 100*(irel * (1+sum(sr)) - sr) # formula 6 in Sugita 2007
    # thats it

  }, simplify=TRUE)

  cm1 <- matrix(unlist(love1),nrow=nrow(loc))
  rownames(cm1) <- loc[[1]]
  lovemean1 <- apply(cm1, 1, mean)
  lovesd1 <- apply(cm1, 1, sd)



  #########################################
  # 2. LOVE run with input from first round
  # for the local composition
  #########################################

  love2 <- replicate(n=repeats,{

    # add random error to PPEs (using mean and stanadard deviation), 'while loop' is to avoid negative values
    PPEs <- rnorm_reveals(length(params[["PPEs"]]), params[["PPEs"]], params[["PPE.errors"]])
    #calculate the s roof from Sugita 2007 LOVE formula 7
    sr2 <- sr2(regsample, lovemean1, params, tBasin = "peatland", dBasin=dBasin, rsap=rsap, dwm=dwm)
    # randomize local counts
    locsampler <- rmultinom(1, sum(locsample), locsample/sum(locsample))

    i <- locsampler/pro
    irel <- i/sum(i)

    fin <- 100*(irel * (1+sum(sr2)) - sr2) # formula 6 in Sugita 2007

    # thats it

  }, simplify=TRUE)


  cm2 <- matrix(unlist(love2),nrow=nrow(loc))
  rownames(cm2) <- loc[[1]]
  lovemean2 <- apply(cm2, 1, mean)
  lovesd2 <- apply(cm2, 1, sd)



  return(as.data.frame(cbind(lovemean1, lovesd1, lovemean1+lovesd1,lovemean2, lovesd2, lovemean2 + lovesd2)))

}


###################################################################
# calculate s roof in Shinyas way
###################################################################

sr <- function(regsample, params, tBasin, dBasin, rsap, dwm){
  # calculates  S_roof from Sugita 2007 LOVE, formula 7
  # for that we need two dispersal factors, one for r to RSAP and one for RSAP to zmax

  # regsample -randomized regional counts
  # params - parameter file
  # tBasins - either "lake" or "peatland"
  # jetzt mal REVEALS laufen lassen..., für ein sample, ohne error
  # aufpassen dass auch in Spalte 1 eine Zahl steht, weil sonst by apply in REVEALS alles zu
  # character wird!

  #dwm <- tolower(dwm) # just avoid confusion
  j <- t(as.matrix(c(1974, regsample))) # just add a year, because REVEALS expects the first number to be some age
  cover <- DISQOVER::REVEALSinR(pollen=j, as.data.frame(params), dwm=dwm, tBasin, dBasin, repeats = 1, verbose = FALSE)

  # hier nun mal provisorisch das Ergebnis - muss sicher auch anders
  vr <- cover[ , grepl("meansim", names(cover) )] / 100  #the v roof needed in Sugitas formula 7

  # produce the local deposition factor (under the Bruch in Sugita 2007 LOVE formel 7)
  f_loc <- do.call(rbind,lapply(params[["fallspeed"]], DispersalFactorK, tBasin, dBasin, dwm, regionCutoff=rsap))

  # produce the regional S_dach (über dem Bruch in Sugita 2007 LOVE formel 7)
  f_reg <- do.call(rbind,lapply(params[["fallspeed"]], DispersalFactorK, tBasin, dBasin=rsap, dwm, regionCutoff = 1e5))

  # und nun alles unterm Bruch
  # das ist totalproportional  contribution distance weighted vegetation to total
  # pollen deposition - under the assumption that composition inside RSAP is the same as outside - Shinyas Magic!
  s_bot <- sum(f_loc*vr) # ub für unterm Bruch

  # und nun den Zähler (regionaler cover * deposition factor)
  s_top <- f_reg * vr
  s <- s_top/s_bot

  # nun ist alles da...

}

###################################################################
# calculate s roof with local composition from first round
###################################################################

sr2 <- function(regsample, lovemean1, params, tBasin, dBasin, rsap, dwm){
  # calculates  S_roof from Sugita 2007 LOVE, formula 7
  # for that we need two dispersal factors, one for r to RSAP and one for RSAP to zmax

  # regsample -randomized regional counts
  # params - parameter file
  # tBasins - either "lake" or "peatland"
  # jetzt mal REVEALS laufen lassen..., für ein sample, ohne error
  # aufpassen dass auch in Spalte 1 eine Zahl steht, weil sonst by apply in REVEALS alles zu
  # character wird!

  #dwm <- tolower(dwm) # just avoid confusion
  j <- t(as.matrix(c(1974, regsample))) # just add a year, muss in REVEALS noch anders werden...
  cover <- DISQOVER::REVEALSinR(pollen=j, as.data.frame(params), dwm=dwm, tBasin, dBasin, repeats = 1, verbose = FALSE)

  # get the mean results
  vr <- cover[ , grepl("meansim", names(cover) )] / 100  #the v roof needed in Sugitas formula 7

  # produce the local deposition factor (under the Bruch in Sugita 2007 LOVE formel 7)
  f_loc <- do.call(rbind,lapply(params[["fallspeed"]], DispersalFactorK, tBasin, dBasin, dwm, regionCutoff=rsap))

  # produce the regional S_dach (über dem Bruch in Sugita 2007 LOVE formel 7)
  f_reg <- do.call(rbind,lapply(params[["fallspeed"]], DispersalFactorK, tBasin, dBasin=rsap, dwm, regionCutoff = 1e5))

  # und nun alles unterm Bruch
  # pollen deposition - now calculated with composition produced in first round
  s_bot <- sum(f_loc*lovemean1) # ub für unterm Bruch

  # und nun den Zähler (regionaler cover * deposition factor)
  s_top <- f_reg * vr
  s <- s_top/s_bot


}
