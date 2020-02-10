#' REVEALSinR
#'
#' The REVEALS model (Regional Estimates of VEgetation Abundance from Large
#' Sites), introduced by Shinya Sugita in 2007, aims to translate pollen data
#' from large sites into regional vegetation composition. The model applies
#' relative pollen productivity estimates (known as PPEs or RPPs) to account for
#' the productivity bias and a "dispersal and deposition factor" K to account
#' for the dispersal bias in pollen data. K is calculated with a dispersal
#' model. REVEALS does not account for the homogeneity bias; we thus have to
#' assume that the vegetation of the pollen source area has been homogeneous.
#' The present implementation of REVEALS in R adds full flexibility with respect
#' to the underlying dispersal model. REVEALSinR by default uses a Lagrangian
#' stochastic (LS) dispersal model (Kuparinen et al. 2007) because this model
#' predicts the pollen dispersal pattern more realistically than the simple
#' Gaussian models used hitherto. Model selection has large effects on the
#' results. REVEALSinR provides further options for testing. To arrive at error
#' estimates, model runs are repeated (by default 1000 times) with random error
#' added in pollen data and PPEs during each model run. To account for pollen
#' deposition in lakes, a lake model is included.
#'
#' REVEALSinR applies calculations on a sequence of samples, such as a pollen
#' record. When needed, the REVEALS function takes calculations for a single
#' sample.
#'
#' @param pollen data frame with raw pollen counts. Counts are arranged as
#'   pollen taxa in columns and samples in rows. The first colum includes sample
#'   ages or depth, i.e. only numeric values! The first row includes taxon
#'   names. Taxon names cannot include spaces! See 'pollenTS' in the 'Tiefer
#'   See' data set to as example. For each pollen taxon, fallspeed and PPEs have
#'   to be provided in params
#' @param params data frame of parameters. Includes 'fallspeed' of pollen (in m
#'   s-1), relative pollen productivity estimates ('PPE') and their standard
#'   error. Pollen taxa are given in rows, parameters in colums.  Column names
#'   must be "species", "fallspeed", "PPE", and "PPE.error". See 'paramsTS' in
#'   the 'Tiefer See' data set as example. For each pollen taxon in the pollen
#'   data, one record of parameters with exactly the same name is required.
#'   Params may include more taxa, also in different order, so that a standard
#'   list may be established.
#' @param dwm distance weighting method. The following methods are implemented:
#'   'lsm unstable', 'gpm unstable', 'gpm neutral'. 'lsm' refers to the the
#'   Lagrangian stochastic model presented by Kuparinen et al. 2007. 'gpm'
#'   refers to the Gaussian plume model, based on Sutton's equations, used by
#'   Prentice (1985; Quaternary Research 23: 76-86). Further methods can be
#'   added.
#' @param tBasin type of basin, either 'lake' oder 'peatland'
#' @param dBasin basin diameter in meter
#' @param n number of model runs per time slice, by default 1000
#' @param regionCutoff	diameter of the reconstruction region in meters, by
#'   default 100000
#' @param ppefun function to randomise these parameters, defaults to
#'   rnorm_reveals
#' @param pollenfun	function to randomise pollen counts, for instance rMaher or
#'   rmultinom_reveals
#' @param writeresults boolean, indicates whether to write results to csv-files
#'   or to return them on the command-line
#' @param verbose	boolean, indicating whether to print progress messages
#'
#' @return REVEALS returns a list with details of the reconstructed vegetation
#'   composition, calculated from the n model runs for each sample
#' @return ...meansim:	mean
#' @return ...mediansim: median
#' @return ...q90sim:	90\% quantile
#' @return ...q10sim: 10\% quantile
#' @return ...sdsim: standard deviation
#'
#' @author Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
#' @author Volkmar Liebscher <volkmar.liebscher@uni-greifswald.de>
#' @references Theuerkauf M., Couwenberg J., Kuparinen A., Liebschar V. (2016) A
#'   matter of dispersal: REVEALSinR introduces state-of-the-art dispersal
#'   models to quantitative vegetation reconstruction. Vegetation History and
#'   Archaeobotany 25: 541-553. https://doi.org/10.1007/s00334-016-0572-0
#'
#' @examples ## Application on a core section from Lake Tiefer See,
#' ## covering the period 1880 t0 2010 CE:
#' data("Tiefer_See", package = "disqover")
#'
#' ## inspect the sample data
#' pollenTS
#' paramsTS
#'
#' ## and apply REVEALS...
#' cover <- REVEALSinR(pollenTS, paramsTS, dwm='lsm unstable', tBasin='lake', dBasin=600)
#'
#' ## plot results as stratigraphic diagram ('pollen diagram')
#' plotREVEALS(cover, type = "single", ytitle = "Age CE")
#'
#' ## plot results as stacked area diagram
#' plotREVEALS(cover, type = "stacked", ytitle = "Age CE")
#'
#' ## load pollen data from an Excel file (e.g.'Pollen_Tiefer_See_1880-2010.xlsx')
#' datafile <- system.file('extdata/reveals', 'Pollen_Tiefer_See_1880-2010.xlsx', package = 'disqover') # with data file from the package
#' datafile <- 'YourFileFolder/Pollen_Tiefer_See_1880-2010.xlsx' # with the file in YourFileFolder
#' pollenTS <- as.data.frame(do.call(read_excel,args=c(list(datafile,sheet = 1)))) # function reads first worksheet
#'
#' ## load parameters from an Excel file (e.g.'pf_PPE_MV2015.xlsx')
#' datafile <- system.file('extdata/reveals', 'pf_PPE_MV2015.xlsx', package = 'disqover') # with data file from the package
#' datafile <- 'YourFileFolder/pf_PPE_MV2015.xlsx' # with the file in YourFileFolder
#' paramsTS <- as.data.frame(do.call(read_excel,args=c(list(datafile,sheet = 1)))) #function reads first worksheet
#'
#' @export


#########################################
##### REVEALSinR - the wrapper function
#########################################

REVEALSinR = function(pollen, params, tBasin, dBasin, dwm = "lsm unstable", n = 1e3, regionCutoff = 1e5,
                      ppefun = rnorm_reveals, pollenfun = rmultinom_reveals,
                      writeresults = FALSE, verbose = TRUE){

  # check main parameters
  # lowercase for better handling...
  dwm <- tolower(dwm); tBasin <- tolower(tBasin)
  names(params) <- tolower(names(params))

  # just in case, make data frame (readxl delivers a different format)
  pollen <- as.data.frame(pollen)
  params <- as.data.frame(params)

  # and change structure of params for better handling
  row.names(params) <- params[,1]
  params <- params[,-1]

  if (!checkREVEALS(pollen = pollen,  params = params, tBasin = tBasin, dBasin = dBasin, dwm = dwm,  n = n,
                    regionCutoff = regionCutoff, writeresults = writeresults, verbose = verbose))
    stop ("Some parameter is not suitable - check and correct")

  # select required species from parameters
  # so there may be more pollen types in params than in the pollen data, and in different order
  params <- params[colnames(pollen[,-1]),]

  # extract PPEs and vg (fall speed)
  ppes <- params[,c('ppe','ppe.error')]
  vg <- params[,'fallspeed']

  # calculate deposition factor K depending on basin size/type and dispersal model
  if(verbose) { message("REVEALS:\tCalculating deposition factor for each species ...")}
  deposition <- do.call(rbind,lapply(vg, DispersalFactorK, tBasin = tBasin, dBasin = dBasin, dwm = dwm,
                                     regionCutoff = regionCutoff))

  # call REVEALS for each time slice in the pollen record
  # die erste Spalte (Age...) sollte besser nicht mit - nur die Pollendaten
  # gibt sonst Probleme wenn das mal keine Zahl ist, weil der Vektor das Format des ersten bekommt!!
  all<-  apply(pollen, 1, REVEALS,  n = n, ppes = ppes, ppefun = ppefun, deposition = deposition,
               pollenfun = pollenfun, writeresults = writeresults, verbose = verbose)

  #making a data.frame for the results statistics
  results_df<-data.frame(dwm, tBasin,pollen[1],t(sapply(all,FUN=function(x)x$meansim)),t(sapply(all,FUN=function(x)x$mediansim)),t(sapply(all,FUN=function(x)x$q90sim)),t(sapply(all,FUN=function(x)x$q10sim)),t(sapply(all,FUN=function(x)x$sdsim)))

  #with the right names
  colnames(results_df)<-c("Distance.weighting", "Basin.type",names(pollen)[1],paste(names(pollen[-1]),rep(c('mean','median','q90','q10','sd'),each=length(pollen[1,])-1),sep="."))

  #writing them if necessary
  if(writeresults)   {
    nn<-names(results_df)
    #quantiles
    k<-c(1:4,grep(x=nn,pattern='mediansim'),grep(x=nn,pattern='q10sim'),grep(x=nn,pattern='q90sim'))
    write.table(results_df[,k],file="ResultStatistics_confidence_limits.csv", append = FALSE, row.names=FALSE, col.names = TRUE, sep=";", dec=",")

    #mean+sd
    k<-c(1:4,grep(x=nn,pattern='meansim'),grep(x=nn,pattern='sdsim'))
    write.table(results_df[,k],file="ResultStatistics_SD.csv",  append = FALSE, row.names=FALSE, col.names = TRUE,  sep=";", dec=",")
  }


  invisible(results_df)
}

#############################################
##### REVEALS - applying the REVEALS model to each time slice
#############################################


REVEALS = function(pollen, n = n, ppes, ppefun = disqover::rnorm_reveals, deposition,
                   pollenfun = disqover::rmultinom_reveals, writeresults = writeresults, verbose = verbose){

  # print time slice to indicate progress
  if (verbose)  {message(paste('REVEALS:\tcomputing for age slice\t',pollen[1]))}

  # start calculation round
  coverM<-replicate(n = n,{

    # add random error to ppes (using mean and stanadard deviation)
    ppes[,3]<-ppefun(length(ppes[,'ppe']), ppes[,'ppe'], ppes[,'ppe.error'])

    # add random error to pollen data (by default: random drawing with rnorm from the given composition)
    pollenPer <- pollenfun(n = 1, pollen = pollen[-1])

    # prepare vector dispersal x PPE
    disp_ppe<-deposition*ppes[,3]

    # calculate cover (cf. REVEALS model formula B1 in Sugita 2007 REVEALS paper)
    s <- sum(pollenPer/disp_ppe)
    round(100*pollenPer/(disp_ppe * s),digits=3)

  },simplify=TRUE)

  # convert results into a matrix
  rownames(coverM) <- rownames(pollen[-1])

  # write complete data file for each time slice
  if(writeresults) {
    fileName <- paste("results_cover", as.character(pollen[1]),".csv",sep="")
    write.csv2(t(coverM), fileName)
  }

  # calculate mean
  coverMean <- apply(coverM, 1, mean)

  # calculate median
  coverMedian <- apply(coverM, 1, median)

  # calculate standard deviation
  coverSD <- apply(coverM, 1, sd)

  # quantiles
  coverq10 <- apply(coverM,1,function(u) quantile(probs=.1,x=u))
  coverq90 <- apply(coverM,1,function(u) quantile(probs=.9,x=u))

  return(invisible(list(meansim=coverMean,sdsim=coverSD,mediansim=coverMedian,q90sim=coverq90,q10sim=coverq10)))
}




###################################################
##### Check parameters
###################################################



checkREVEALS <- function(pollen,  params, tBasin, dBasin, dwm,  n, regionCutoff, writeresults, verbose){

  # neu...
  if ((dBasin/2) > regionCutoff) stop("Stop: Basin larger than the size of the region (regionCutoff) - that does'nt work!")

  # test whether all pollen taxa are in the parameter list (their order is not important)
  if(!all(colnames(pollen[,-1]) %in% rownames(params), na.rm = FALSE))
    stop("STOP: Not all pollen taxa found in the parameter list")

  # Ages/depth numeric?
  if(!is.numeric(pollenTS[,1]))
    stop("STOP: some or all values in first column of pollen record (ages/depths) not numeric")

  # check distance weighting method
  if (!isTRUE((dwm=="lsm unstable") | (dwm=="gpm neutral") | (dwm=="gpm unstable") | (dwm=="1overd")))
    stop("distance weighting method not defined; should be 'LSM unstable', 'GPM neutral', 'GPM unstable' or '1oved' ")

  # check basin type
  if (!isTRUE((tBasin=="peatland") | (tBasin=="lake")))
    stop("basin type ('tBasin') not defined; should be 'peatland' or 'lake'")

  # check basin diameter
  if (!isTRUE(dBasin == floor(dBasin))) stop("basin size ('dBasin') should be an integer value")
  if (dBasin < 10) stop("basin diameter ('dBasin') should be at least 10 m")

  # check number if iterations (should be larger than 1000)
  if (!isTRUE(n == floor(n))) stop("'n' should be an integer value")
  if (n < 1000) message("Warning: for sensible error estimates, 'n' should be at least 1000") # kann ruhig auch mal kleiner sein...

  #fallspeed column present?
  if(!'fallspeed' %in% tolower(names(params)))
    stop("STOP: no 'fallspeed' column in the parameters ")

  # values between 0.01 and 0.15?
  if(max(params$fallspeed)>0.15) stop ("fallspeed(s) too high (>0.15 m/s), please check")
  if(min(params$fallspeed)<0.01) stop ("fallspeed(s) too low (<0.01 m/s), please check")

  # ppe column present?
  if(!'ppe' %in% tolower(names(params)))
    stop("STOP: no 'ppe' column in the parameters ")

  # PPE.errors column present?
  if(!'ppe.error' %in% tolower(names(params)))
    stop("STOP: no 'ppe.error' column in the parameters ")

  return(TRUE)
}


#########################################
##### add error in pollen data
#########################################

# default: add errors to pollen data through random drawing
rmultinom_reveals<-function(n=1,pollen,...)
{
  100*rmultinom(n, sum(pollen), pollen/sum(pollen))/sum(pollen, na.rm = TRUE)
}


# alternative: add error to pollen data following Maher (1972)
rMaher <- function(n=1,pollen,...){

  countSum<-sum(pollen[-1])
  pDat <- pollen[-1]/countSum

  ranVec<-rnorm(length(pollen)-1, mean=0, sd=1) #random Vector with mean of zero and 1 standard deviation
  i1 <- pDat+((ranVec^2)/(2*countSum))
  i2 <- (pDat*(1-pDat)/countSum)+((ranVec^2)/(4*countSum^2))
  denom <- 1+((ranVec^2)/countSum)
  ranPol <- 100*(i1+(ranVec*sqrt(i2)))/denom
  return(ranPol)
}


#########################################
##### add error in ppes
#########################################

# changing the PPE's randomly
# 'while loop' is to avoid negative values
rnorm_reveals <- function(n , mean=rep(1, n), sds = rep(1 ,n))
{x <- 0 * mean - 1
 #iterate till positivity is established
 while(any(x < 0)){x <- rnorm(n = n, mean = mean, sd = sds)}
 return(x)
}





