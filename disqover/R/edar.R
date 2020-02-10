#' Extended Downscaling Approach (EDA)
#'
#' The Extended Downscaling Approach (EDA, Theuerkauf and Couwenberg 2017) aims
#' to detect past vegetation patterns and communities within landscapes. The EDA
#' uses iterative forward modelling to fit vegetation composition to robust
#' landscape patterns by comparing simulated with observed pollen deposition.
#' The approach employs a set of pollen records, preferably from medium sized to
#' large lakes or peatlands, as well as environmental data extracted from e.g.
#' soil maps or digital elevation models. The EDA calculations need time, so
#' plan several hours.
#'
#' @param edav a list with pollen data, environmental data and parameters, best
#'   produced with the import function edaimportxl()
#' @param dwm distance weighting method, defaults to 'lsm_unstable' (Lagrangian
#'   stochastic model), further models are the Gaussian plume model ('gpm
#'   neutral' and 'gpm unstable') and the non-parametric model '1overd'.
#' @param n number of EDA rounds. To produce error estimates, EDA calculations
#'   are repeated 'n' times with noise added in pollen counts and/or PPEs. Error
#'   estimates requires a minimum of 51 rounds. The minimum of 'n' is 2
#' @param pollenError adds error in pollen data through resampling pollen
#'   counts, defaults to TRUE
#' @param ppeError adds error in PPEs on the basis of the given standard error.
#'   Defaults to TRUE
#' @param itermax DEoptim control parameter: The maximum iteration (population
#'   generation) allowed. Defaults to 20000
#' @param reltol  DEoptim control parameter: relative convergence tolerance. The
#'   algorithm stops if it is unable to reduce the value by a factor of reltol *
#'   (abs(val) + reltol) after steptol steps. Here defaults to 0.01
#' @param steptol  DEoptim control parameter: see reltol. Here defaults to 500
#' @param trace  DEoptim control parameter: Positive integer or logical value
#'   indicating whether printing of progress occurs at each iteration. The
#'   default value is TRUE. If a positive integer is specified, printing occurs
#'   every trace iterations
#' @param parallelType  DEoptim control parameter: Defines the type of
#'   parallelization to employ, if any. 0: This uses DEoptim on one only one
#'   core. 1: The default, this uses all available cores, via the parallel
#'   package, to run DEoptim. If parallelType=1, then the packages argument and
#'   the parVar argument need to specify the packages required by the objective
#'   function and the variables required in the environment, respectively. 2:
#'   This uses the foreach package for parallelism; see the sandbox directory in
#'   the source code for examples. If parallelType=2, then the foreachArgs
#'   argument can pass the options to be called with foreach
#'
#' @return The reconstructed local vegetation composition
#' @return ...median: median cover of each taxon on each environmental class
#' @return ...mean: mean cover of each taxon on each environmental class
#' @return ...quantile10: 10\% quantile of the cover of each taxon on each
#'   environmental class
#' @return ...quantile90: 90\% quantile of the cover of each taxon on each
#'   environmental class
#' @return ...bestmem: optimized cover values from each run (n = number of taxa
#'   x number of environmental classes)
#' @return ...bestval: minimum of the target value in each run
#' @return ...iter: number of iterations in each run
#'
#' @author Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
#' @references Theuerkauf M. and Couwenberg J. (2017) The extended downscaling
#'   approach: A new R-tool for pollen-based reconstruction of vegetation
#'   patterns. The Holocene 27(8): 1252–1258.
#'   https://doi.org/10.1177/0959683616683256
#'
#' @examples # Application with example data from Theuerkauf & Couwenberg (2017)
#' # load data set "Scenario C"
#' data(eda_tc17c, package = "disqover")
#'
#' # view structure and elements of the data set
#' names(eda_tc17c)
#' names(eda_tc17c$env)
#' str(eda_tc17c)
#'
#' # and run EDA with standard parameters
#' e <- edar(edav = eda_tc17c)
#'
#' # edaimportxl() helps to prepare the data from Excel files. It needs three files:
#' # pollen data, environmental data and taxon parameters. The following example shows the principle
#' # with external data from the package. First define file names:
#' pollenfile <- system.file("extdata/eda", "eda_tc17c_pollen.xlsx", package = "disqover")
#' envfile <- system.file("extdata/eda", "eda_tc17c_environment.xlsx", package = "disqover")
#' paramsfile <- system.file("extdata/eda", "eda_tc17c_params.xlsx", package = "disqover")
#'
#' # now use edaimportxl()
#' eda_tc17c <- edaimportxl(pollenfile, envfile, paramsfile)
#'
#' @export


#########################################
##### edar - the wrapper function
#########################################


edar = function(edav, dwm = "LSM unstable", n = 100, pollenError = TRUE, ppeError = TRUE,
                itermax = 20000, reltol = 0.01, steptol = 500, trace = FALSE, parallelType = 1, verbose = TRUE, u = 3){

  # requires compiler
  if (!requireNamespace("compiler", quietly = TRUE)) {
    stop("Package compiler needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # requires DEoptim
  if (!requireNamespace("DEoptim", quietly = TRUE)) {
    stop("Package DEoptim needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # compile asa.f to be faster
  asa.f.c <- compiler::cmpfun(asa.f)

  # for better handling...
  dwm <- tolower(dwm)
  colnames(edav$params) <- tolower(colnames(edav$params))

  # check input
  checkeda(edav = edav, dwm = dwm, dBasin = dBasin, n = n, pollenError = pollenError,
           ppeError = ppeError, itermax = itermax, reltol = reltol, steptol = steptol,
           trace = trace, parallelType = parallelType, verbose = verbose, u = u)


  # select required species from parameters
  edav$params <- edav$params[row.names(edav$pollen),]

  # pollen <- as.matrix(pollen) # needed?

  dwm<-tolower(dwm) # all in lowercase to make the test easier

  # extract fall speed
  vg <- edav$params$fallspeed

  # divide cover in each ring by ring area
  radii <- as.double(c(0, edav$radii)) # extract radii, add zero for first ring, unname would remove names

  areaCircles <- pi * (radii^2); areaRings <- abs(diff(areaCircles))  # calculate ring areas
  envr <- lapply(edav$env, function(x) x/areaRings) # geht - aber auch richtig?

  # produce distance weighted cover of each envionmental unit/class
  # for each taxon (=fall speed) seperately, i.e. each list element corresponds to one taxon
  envdw <- lapply(vg, dwpa, dwm = dwm, envdata = envr, radii = radii, u = u)
  message("distance weighting completed, optimization starts, now be patient...")

  # extract some parameters
  nTaxa <- length(edav$params[[2]])  # read number of plant taxa
  nClasses <- length(edav$env)  # read number of landcover classes
  nLakes <- ncol(edav$pollen)

  # prepare status message
  message('Run Nr:')

  # set bounaries for optimization
  lower <- rep(0, nTaxa*nClasses)
  upper <- rep(1, nTaxa*nClasses)

  ##################################################################
  ### eda core function
  ##################################################################

  eda <- function(){

    # to deal with NAs, set them zero but remember
    nas <- is.na(edav$pollen)
    edav$pollen[nas] <- 0

    if (isTRUE(pollenError)) {
      b <- apply(edav$pollen, 2, multidraw)
      pollenPer <- 100*prop.table(b, 2)
    } else {
      pollenPer <- 100*prop.table(as.matrix(edav$pollen),2)
    }

    # add the NAs again
    pollenPer[nas] <- NA

    # produce PPEs with/without random error
    # refernce Taxon wird jetzt noch mit verändert!! Mist
    # kann jetzt auch noch negativ werden - sollte nicht
    if (isTRUE(ppeError)) {
      ppe<-rnorm(length(edav$params[,c('ppe')]), edav$params[,c('ppe')], edav$params[,c('ppe.error')])
      while (min(ppe) <= 0) {
        ppe<-rnorm(length(edav$params[,c('ppe')]), edav$params[,c('ppe')], edav$params[,c('ppe.error')])
        }
    } else {
      ppe<-edav$params[,c('ppe')]
    }

    oX <- DEoptim::DEoptim(asa.f.c, lower, upper, ppe = ppe, envdw = envdw, pollenPer = pollenPer, nTaxa = nTaxa,
                          DEoptim::DEoptim.control(NP = 11*nTaxa*nClasses,
                          strategy = 2, itermax = itermax, reltol = reltol, steptol = steptol, trace = trace,
                          parallelType = parallelType, VTR = 0))

    # the future:
    # oX <- RcppDE::DEoptim(asa.f.c, lower, upper, ppe = ppe, envdw = envdw, pollenPer = pollenPer, nTaxa = nTaxa,
    #      RcppDE::DEoptim.control(VTR = 0, strategy = 2, NP = 11*nTaxa*nClasses, itermax = itermax, trace = trace))

    return(oX)
  }


  # prepare result variables for bestmem values
  bm <- matrix(nrow = nTaxa*nClasses, ncol = n)
  bma <- vector("list", nClasses)
  bmp <- vector("list", nClasses)

  cmedian <- matrix(nrow = nTaxa, ncol = nClasses)
  cmean <- matrix(nrow = nTaxa, ncol = nClasses)
  cq10 <- matrix(nrow = nTaxa, ncol = nClasses)
  cq90 <- matrix(nrow = nTaxa, ncol = nClasses)

  bv <- matrix(1, ncol = n)
  it <- matrix(1, ncol = n)

  for (i in seq_along(1:n)){
    if (verbose) message(i)
    res <- eda()
    bm[,i] <-unlist(res$optim$bestmem)
    bv[,i] <-unlist(res$optim$bestval)
    it[,i] <-unlist(res$optim$iter)
  }

  # summarize the results
  for (i in seq_along(1:nClasses)){
    f <- ((i-1)*nTaxa+1); t <- (i*nTaxa)
    bma[[i]] <- bm[f:t, ]
    bmp[[i]] <- 100*prop.table(bma[[i]], 2)
    cmedian[ ,i] <- (apply(bmp[[i]], 1, median))
    cmean[ ,i] <- (apply(bmp[[i]], 1, mean))

    # quantiles only with at least 51 runs
    if (n>50){
      cq10[ ,i] <- (apply(bmp[[i]], 1, function(u) quantile(probs = .1, x = u)))
      cq90[ ,i] <- (apply(bmp[[i]], 1, function(u) quantile(probs = .9, x = u)))
    }
  }

  if (n<=50){
    colnames(cmedian)<-names(edav$env)
    rownames(cmedian)<-rownames(edav$pollen)

    colnames(cmean)<-names(edav$env)
    rownames(cmean)<-rownames(edav$pollen)

    rl <- list(round(cmedian,3), round(cmean,3), bm, bv, it)
    names(rl) <- c("median", "mean", "bestmem", "bestval", "iter")
  }

  # 10% and 90% quantiles only when sufficient run (50)
  if (n>50){
    colnames(cmedian)<-names(edav$env)
    rownames(cmedian)<-rownames(edav$pollen)

    colnames(cq10)<-names(edav$env)
    rownames(cq10)<-rownames(edav$pollen)

    colnames(cq90)<-names(edav$env)
    rownames(cq90)<-rownames(edav$pollen)

    rl <- list(round(cmedian,3), round(cmean,3), round(cq10,3), round(cq90,3), bm, bv, it)
    names(rl) <- c("median", "mean","quantile10", "quantile90", "bestmem", "bestval", "iter")
  }

  return(rl)
}



##################################################################################
# produde pollen percentages with (using multinomial drawing) or without random error
##################################################################################

multidraw <- function(a){

  rmultinom(1, sum(a), a/sum(a))
}


##################################################################################
##### asa.f   calculates distance between modeled and empiric pollen percentages
##################################################################################

# parameters:
# cover: a vector that gives the abundance of each plant taxon on each soil/substrate types --> to be optimized
# ppe: a vector with  pollen productivity of each pollen taxon
# envdw - a list of matrics where each matrix gives the distance weighted cover
# of one soil/substrate type in one (of many) rings around each lake
# pollenPer: observed pollen percentages for a number of taxa (=nTaxa) in a number of lakes
# cover <- upper
asa.f = function(cover, ppe, envdw, pollenPer, nTaxa){

  cp <- matrix(cover*ppe, nrow = nTaxa) # cover times PPE as matrix

  # multiply plant proportions (to be optimized) and PPE with distance weighted land cover
  pollenInflux <- list()
  envdwcp <- envdw
  for (i in 1:nTaxa){
    envdwcp [[i]] <-  cp[i,]*envdw[[i]]
    pollenInflux[[i]] <- colSums(envdwcp [[i]])
  }

  influx <- do.call(rbind, pollenInflux) # ist das schnell?

  # calculate modeled pollen percentages (with NAs allowed)
  pollenMod <- 100*prop.table(influx, 2)

  # to consider NAs...
  nas <- is.na(pollenPer)
  pollenPer[nas] <- 0

  dif <- ((pollenMod-pollenPer)^2)/((pollenPer+1)*ppe)
  dif[nas] <- NA
  asa <- sum(dif, na.rm = TRUE)

  # calculate distance between modelled and empiric pollen data
  # asa <- sum(((pollenMod-pollenPer)^2)/((pollenPer+1)*ppe))

  return(asa)
}


##################################################################################
##### dwpa   produces a distance weighted matrix of land cover for each taxon
##################################################################################

# parameters:
# vg - fall speed of pollen
# dwm - distance weighting method
# envdata - the environmental data
# radii - a list of radii
# u - wind speed (only needed for the gpm models, by default = 3)

dwpa= function(vg, dwm, envdata, radii, u){

  airborne <-  DispersalFactorK(vg=vg, tBasin="peatland", dBasin = 2*radii, dwm = dwm, regionCutoff = 1e7, u)
  influxRing <- abs(diff(airborne))

  lcSum <- list()
  for (i in 1:length(envdata)){
    # multiply landcover with distance weightin factor --> weighted land cover
    envdata[[i]] = envdata[[i]] * influxRing

    # row sums --> sum up landcover in all rings around a lake
    lcSum[[i]] <- apply(envdata[[i]], 2, sum)
  }

  lcSum <- do.call(rbind, lcSum)
  return(lcSum)
}


###################################################
##### Check parameter file
###################################################

checkeda <- function(edav = edav, dwm = dwm, dBasin = dBasin, n = n, pollenError = pollenError,
                     ppeError = ppeError, itermax = itermax, reltol = reltol, steptol = steptol,
                     trace = trace, parallelType = parallelType, verbose = verbose, u = u){

  ### first check edav with the pollen, parameters and environmental data

  # test whether all pollen taxa are in the parameter list (their order is not important)
  if(!all(rownames(edav$pollen) %in% rownames(edav$params), na.rm = FALSE))
    stop("STOP: Not all pollen taxa are in the parameter list")

  # test whether the columns 'fallspeed', 'ppe', and 'ppe.error' are included in the parameter file
  if(!'fallspeed' %in% tolower(colnames(edav$params)))
    stop("STOP: no 'fallspeed' column in the parameters ")

  if(!'ppe' %in% tolower(colnames(edav$params)))
    stop("STOP: no 'ppe' column in the parameters ")

  if(!'ppe.error' %in% tolower(colnames(edav$params)))
    stop("STOP: no 'ppe.error' column in the parameters ")

  # check range of fallspeeds (limited through LS model, may be extended upon request)
  if(max(edav$params$fallspeed)>0.15) stop ("one or more fallspeed values too high (>0.15 m/s), please check")
  if(min(edav$params$fallspeed)<0.01) stop ("one or more fallspeed values too low (<0.01 m/s), please check")

  # check dimension in environmental data
  r <- lapply(edav$env, nrow)
  if (!all(r == r[[1]])) stop("Warning: number of radii (=rows) in environmental data not equal")
  c <- lapply(edav$env, ncol)
  if (!all(c == c[[1]])) stop("Warning: number of sites (=colums) in environmental data not equal")


  ### and now the parameters

  # dwm - distance weighting method
  if (!isTRUE((dwm=="lsm unstable") | (dwm=="gpm neutral") | (dwm=="gpm unstable"))) stop("distance weighting method not defined; should be 'LSM unstable', 'GPM neutral' or 'GPM unstable'")

  # n - calculation rounds
  if (!isTRUE(n == floor(n))) stop("'n' should be an integer value")
  if (n < 2) stop("'n' (number of rounds) should be at least 2")
  if (n < 51) message("Warning: n < 51. No error estimates  will be produced with less than 51 rounds")

  # pollenError
  if (!is.logical(pollenError)) stop("'pollenError' should be TRUE or FALSE")

  # ppeError
  if (!is.logical(ppeError)) stop("'ppeError' should be TRUE or FALSE")

  # itermax
  if (!isTRUE(itermax == floor(itermax))) stop("'itermax' should be an integer value")
  if (itermax < 2000) message("Warning: 'itermax' rather small, use higher value (e.g. 20000) for optimal results")

  # reltol
  # ???

  # steptol
  # ???

  # trace
  if (!is.logical(trace)) stop("'trace' should be TRUE or FALSE")

  # parallelType - not used in RcppDE
  if (!(parallelType %in% c(0,1,2))) stop("'parallelType' out of range (0, 1 or 2)")

  # verbose
  if (!is.logical(verbose)) stop("'verbose' should be TRUE or FALSE")

  # u (wind speed)
  if (!is.numeric(u)) stop("'u' (windspeed) should be numeric")
  if (u < 0) stop("'u' (windspeed) cannot be negative")
  if (u > 20) message("'u' (windspeed) seems overly high")

  return(TRUE)
}




###################################################
##### import functions: Excel
###################################################



#' Preparing data for the Extended Downscaling Approach with edar()
#'
#' Applying the Extended Downscaling Approach (EDA) with the edar() function
#' requires a number of data: pollen data, environmental data and taxon
#' parameters. All the data has to be provided in a single list. edaimportxl
#' builds this list from three data files. Two of them, the file for pollen data
#' and the file for taxon parameters shall only include one sheet. The file for
#' environmental data includes a number of sheets, one for each class/unit (e.g.
#' soil type). For the structure of each file, see the examples from Theuerkauf
#' & Couwenberg (2017), which are provided as external data.
#'
#' @param pollenfile path and name of the pollen file. This file includes one
#'   sheet with raw pollen counts for each taxon and site (usuall lake). Each
#'   row represents one pollen taxon, each column one site. The first row
#'   includes site names, the first column taxon names. See
#'   'eda_tc17c_pollen.xlsx' as example.
#'
#' @param envfile path and name of the environmental data file. This file
#'   includes spatial information, i.e. the cover of each environmental
#'   unit/class (e.g. each soil type)  within a specific distance from each
#'   site. The area is given as the total area inside the respective full circle
#'   around each site (not in donuts). The file includes one sheet for each
#'   environmental unit. In each sheet, each row represents one radius (in
#'   increasing order), each column one site. The first row includes site names,
#'   the first column radii (as integer or real number). Sheet names should
#'   represent the environmental unit/class. See 'eda_tc17c_environment.xlsx' as
#'   example.
#'
#' @param paramsfile path and name of the parameter file. This file includes
#'   fallspeed of pollen (in m/s), pollen productivity estimates and standard
#'   error of the PPEs. PPEs are region specific corrections factors, so values
#'   from a suitable calibration study are required. PPE calibration is highly
#'   sensitive to dispersal model selection. Therefore, the same underlying
#'   dispersal model has to be applied in PPE calibration and EDA application.
#'   The parameters have to be provided in three columns: 'fallspeed', 'ppe'
#'   and 'ppe.error'. Names of the pollen taxa are represented in row names. See
#'   'eda_tc17c_params.xlsx' as an example. The parameter file may include more
#'   pollen taxa than the pollen data, and in a different order.
#'
#'
#' @return ...edav: a list with input data for edar()
#'
#' @author Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
#' @references Theuerkauf M. and Couwenberg J. (2017) The extended downscaling
#'   approach: A new R-tool for pollen-based reconstruction of vegetation
#'   patterns. The Holocene 27(8): 1252–1258.
#'   https://doi.org/10.1177/0959683616683256
#'
#' @examples # edaimportxl() prepares data for edar()
#' # it requires three files: pollen data, environmental data and taxon parameters
#' # the following example uses external data included in the package
#'
#' # first set path and filenames for the data files
#' pollenfilename <- system.file("extdata/eda", "eda_tc17c_pollen.xlsx", package = "disqover")
#' envfilename <- system.file("extdata/eda", "eda_tc17c_environment.xlsx", package = "disqover")
#' paramsfilename <- system.file("extdata/eda", "eda_tc17c_params.xlsx", package = "disqover")
#'
#' # now use edaimportxl()
#' eda_tc17c <- edaimportxl(pollenfilename, envfilename, paramsfilename)
#'
#' @export

edaimportxl <- function(pollenfile, envfile, paramsfile){
  edav <- list()

  # requires readxl
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package readxl needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # pollen data
  # may contain NA values
  edav$pollen <- readxl::read_excel(path=pollenfile, na = "NA")
  n <- unlist(edav$pollen[1])
  edav$pollen <- as.data.frame(edav$pollen[-1])
  row.names(edav$pollen) <- n

  # parameters
  edav$params <- readxl::read_excel(path=paramsfile)
  n <- unlist(edav$params[1])
  edav$params <- as.data.frame(edav$params[-1])
  row.names(edav$params) <- n

  # environmental data - muss auch noch umegbaut werden, 1. Spalte als Namen
  edav$env <- lapply(readxl::excel_sheets(envfile), readxl::read_excel, path = envfile)
  edav$radii <- unname(unlist(edav$env[[1]][1]))
  edav$env <- lapply(edav$env, function(x) as.data.frame(x[-1])) # convert to data frame
  names(edav$env) <- readxl::excel_sheets(envfile)


  return(edav)
}
