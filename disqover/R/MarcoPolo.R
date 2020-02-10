#' MAnipulating pollen sums to ReCOnstruct POllen of Local Origin
#'
#' MARCO POLO (MAnipulating pollen sums to ReCOnstruct POllen of Local Origin) is a quantitative approach
#' to stand-scale palynology, which translates pollen deposition from very small mires or ponds into
#' local vegetation composition. To that end, MARCO POLO removes those pollen types from the pollen sum
#' whose values are significantly higher than in a neighbouring large basin. The resulting regional pollen sum
#' is free of the disturbing factor of (extra-)local pollen. Based on this sum, comparison with the
#' pollen record from the large basin allows calculating sharp (extra-)local signals.
#' Treating the (extra-)local pollen portion with representation factors (R-values) then produces a
#' quantitative reconstruction of the stand-scale vegetation composition.
#'
#' @param local local pollen data
#' @param regional regional pollen data
#' @param params table of R-values
#' @param assignment table naming the regional record for each local record
#' @return The reconstructed local vegetation composition
#' @author Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
#' @author Volkmar Liebscher <volkmar.liebscher@uni-greifswald.de>
#' @author Almut Mrotzek <almut.mrotzek@gmx.de>
#' @references Mrotzek A., Couwenberg J., Theuerkauf M. and Joosten H. (2017) MARCO POLO –
#'   A new and simple tool for pollen-based stand-scale vegetation reconstruction. The Holocene 27: 321-330.
#'   DOI:10.1177/0959683616660171
#'
#' @examples
#'   #### apply Marco Polo with the sample data from Mrotzek et al. (2017)
#'   data("MP_example", package = "disqover")
#'   a <- MarcoPolo(localMP, regionalMP, paramsMP, assignmentMP)
#'
#'   #### data may be read from seperate files or from a single file with 4 sheets, as in the following example.
#'   #### for data structure see example data (e.g. MarcoPoloExampleData.xlsx)
#'   datafile <- system.file("extdata/marcopolo", "MarcoPoloExampleData.xlsx", package = "disqover")
#'   localMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="local")))
#'   regionalMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="regional")))
#'   paramsMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="params")))
#'   assignmentMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="assignment")))
#'   a <- MarcoPolo(localMP, regionalMP, paramsMP, assignmentMP)
#'
#'   ### plot the results for surface samples
#'   ### prepare the plot
#'   p <- plotMP(a, ytitle = "Site")
#'
#'   ### and write plot to a pdf file in the working directory
#'   pdf("MarcoPolo plot.pdf", height = 12, width = 12)
#'   p
#'   dev.off()
#'
#'
#' @export


#################################################################
### Marco Polo main function
#################################################################

MarcoPolo <- function(local, regional, params, assignment) {

  # checks müssen überarbeitet werden, passen so nicht mehr
  # if  (!(identical(dim(local), dim(regional)))) {
  #   stop(" Interrupted because local and regional have different dimensions...")
  # }
  #
  # # test dimension of params, should be similar to local/regional data set
  # if  (!(nrow(params)+1 == nrow(local))) {
  #   stop(" Interrupted because list of params does not match local/regional data...")
  # }

  # auch checken dass keine Leerzeichen in regional sites, weil die löscht R und dann wird es komisch
  # was checken: local, regional und R-values haben gleiche Länge und gleiche Taxa
  # Anzahl local sites stimmt mit assignments zusammen

  # create data frame
  local.m<-as.matrix(local[,-1]) # local data as matrix
  class(local.m) <- "numeric"
  regional.m<- as.matrix(regional[,-1]) # regional data as matrix
  class(regional.m) <- "numeric"
  v.check <- rep(1,nrow(params)) #vector needed
  pLocal <- rep(0,nrow(params)) # vector needed
  params <- as.data.frame(params) #because read_excel delivers a tribble that behaves differently...

  # start calculation round
  # Tipps for the loop: http://stackoverflow.com/questions/14490182/repeating-a-user-defined-function-using-replicate-or-sapply
  a <- lapply(1:ncol(local.m), function(i) MarcoPolo.main(i, local.m, regional.m, assignment, v.check, pLocal))
  b <- as.data.frame(a)
  c <- b * (params[,2]) #unlist b
  d <- apply(c, 2, function(x) 100*x/sum(x))
  colnames(d) <-colnames(local.m)
  rownames(d) <- (params[,1])

  return(d)
}


#################################################################
### Marco Polo application for 1 local and 1 regional data set
#################################################################

MarcoPolo.main <- function(n, local.m, regional.m, assignment, v.check, pLocal, minSum=20){

  # select name of the respective regional record
  r <- unlist(assignment[n,2])

  # create data frame with 1 sample from local data and respective regional record
  df <- data.frame(local.m[,n], regional.m[,r], v.check, pLocal)

  # provide column names
  names(df)<- c("local", "regional", "check", "pLocal")

  # validate difference between local and regional data
  MarcoPolo.checkDifference <- function(){

    # add Mosimann error to regional pollen values
    reg <- Mosiman(df$regional*df$check, upper=TRUE)

    # add Mosimann error to regional pollen values
    loc <- Mosiman(df$local*df$check, upper=FALSE)

    # count the cases with local values larger than regional values
    nOut <- sum((reg-loc) < 0  )

    # temporary vector of pollen types that remain in the sum
    atest <- df$check*((sign(sign(reg-loc)+0.5)+1)/2)

    # update vector of types used to calculate the pollen sum long as pollen sum remains larger than threshold (by default 20)
    # interupt when new sum smaller than threshold
    if (sum(df$local*atest)>minSum){
      df$check <<- atest
    } else nOut=0

    # interrupt if either no new local types have been detected or the remaining pollen sum is below the threshold
    return(nOut)
  }


  # iteratively identify pollen types with local over-representation
  while (MarcoPolo.checkDifference()>0) {}

  # calculate pollen percentages based on a sum exclduing taxa with local presence
  df$regional <- 100*pPercent(df$regional, df$check)
  df$local <- 100*pPercent(df$local, df$check)

  # calculate "cover" of the local taxa
  df$pLocal <- (df$local-df$regional)*(df$check-1)^2
  df$pLocal <- 100*pPercent(df$pLocal)

  return(df$pLocal)
}




#########################################
##### add Mosiman error range in pollen data
#########################################

# adds error ranges in pollen data following Mosimann (1965) and Maher (1972)
# by default, error range is 1 standard deviation

Mosiman <- function(pollen, upper, SD=1){

  N <- sum(pollen)
  P <- pollen/N
  F <- P+((SD*SD)/(2*N))
  S <- SD * sqrt((P*(1-P)/N)+(SD*SD/(4*N*N)))
  D <- 1+(SD*SD/N)

  if (upper){ result <- 100*(F+S)/D } # Mosimans "H" value

  else { result <- 100*(F-S)/D } # Mosimans "L" value

  return(result)
}



#############################################
### function to calculate pollen percentage
### allows selection
#############################################

pPercent <- function(vec, selection=1){
  vec <- vec/sum(vec*selection)
  return(vec)
}
