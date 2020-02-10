#' plotREVEALS
#'
#' Plots landcover reconstructed with REVEALSinR.
#'
#' @param d the output of REVEALSinR
#' @param type plot type: 'single' for a stratigraphic diagram with seperate
#' columns. Error plotted as 10\% and 90\% qantiles; 'stacked' for stacked area plot.
#' Colours by default taken from the terrain.colors palette
#' @param ytitle  y-axis title, defaults to "Age cal. BP"
#'
#' @return A plot object
#'
#' @author Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
#' @author Triin Reitalu <triinreitalu@gmail.com>
#'
#' @examples ## Apply REVEALSinR and plotting on a core section from Lake Tiefer See:
#'
#' # load 'Tiefer See' data set and apply REVEALSinR
#' data("Tiefer_See")
#' cover <- REVEALSinR(pollenTS, paramsTS, tBasin = "lake", dBasin = 600)
#'
#' # produce stratigraphic plot (like a pollen diagram)
#' # by default, the y-axis is reversed to suit BP dates
#' p <- plotREVEALS(d = output, type = "single", ytitle = "Age CE", reverse = FALSE)
#'
#' # and write to pdf file
#' pdf("reveals_single.pdf", height = 8, width = 12)
#' p
#' dev.off()
#'
#' # now produce stacked area plot
#' p <- plotREVEALS(d = output, type = "stacked", ytitle = "Age CE")
#'
#' # and write to pdf file
#' pdf("reveals_stacked.pdf", height = 12, width = 12)
#' p
#' dev.off()
#'
#' @export


plotREVEALS <- function(d, type = "single", ytitle = "Age cal. BP", reverse = TRUE, error = "quantiles"){

  # Richtung y- axis ändern - je nach Age Spalte (also CE oder BP), oder einfach reverse anbieten

  # plotting a stratigraphic diagramm with seperate colums for each taxon
  # inspired by code from Olivier Blarquez
  # from http://paleoecologie.umontreal.ca/684-2/

  # requires ggplot2 - is now in depends
  # if (!requireNamespace("ggplot2", quietly = TRUE)) {
  #   stop("Package ggplot2 needed for this function to work. Please install it.",
  #        call. = FALSE)
  # }

  # requires forcats
  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop("Package forcats needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #extracting the data
  age <- d[[3]]

  d.median <- d[grep('median', names(d), value=TRUE)]
  names(d.median) <- sub(".median", "", names(d.median))

  d.mean <- d[grep('mean', names(d), value=TRUE)]
  names(d.mean) <- sub(".mean", "", names(d.mean))

  d.10 <- d[grep('q10', names(d), value=TRUE)]
  names(d.10) <- sub(".q10", "", names(d.10))

  d.90 <- d[grep('q90', names(d), value=TRUE)]
  names(d.90) <- sub(".q90", "", names(d.90))

  d.sd <- d[grep('sd', names(d), value=TRUE)]
  names(d.90) <- sub(".sd", "", names(d.sd))

  # now create data frames for median, 10% and 90% quantile
  # warum sortiert wer hier??
  df = data.frame(yr=rep(age,ncol(d.median)),
                per=as.vector(as.matrix(d.median)),
                taxa=as.factor(rep(colnames(d.median),each=nrow(d.median))))

  df.mean = data.frame(yr=rep(age,ncol(d.mean)),
                per=as.vector(as.matrix(d.mean)),
                taxa=as.factor(rep(colnames(d.mean),each=nrow(d.mean))))

  df.10 = data.frame(yr=rep(age,ncol(d.10)),
                   per=as.vector(as.matrix(d.10)),
                   taxa=as.factor(rep(colnames(d.10),each=nrow(d.10))))

  df.90 = data.frame(yr=rep(age,ncol(d.90)),
                   per=as.vector(as.matrix(d.90)),
                   taxa=as.factor(rep(colnames(d.90),each=nrow(d.90))))

  df.sd = data.frame(yr=rep(age,ncol(d.sd)),
                     per=as.vector(as.matrix(d.sd)),
                     taxa=as.factor(rep(colnames(d.sd),each=nrow(d.sd))))

    # add 10% and 90% quantila to main data frame
  if (error == "quantiles"){
    df$mean = df$per
    df$low = df.10$per
    df$up = df.90$per
  } else if (error == "sd"){
    df$mean = df.mean$per
    df$low = df.mean$per - df.sd$per
    df$up = df.mean$per + df.sd$per
  } # sonst Fehlermeldung

  # define a theme without gridlines, etc. and elements that will decrease readability of the diagram
  theme_new <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grids
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size = 10, angle = 90, vjust = 0, hjust = 0.2), # Taxa names
                     strip.background = element_blank(),
                     strip.text.y = element_text(angle = 0),
                     legend.position = "right", panel.border = element_blank(),
                     axis.text.x = element_text(angle=0,hjust=0.5)) # Axis tick label angle

  # and finally the diagrams
  # single taxa
  if (type == "single") {

    p <- ggplot(df) +
      geom_area(aes(yr,per), fill = "grey90") +
      geom_ribbon(aes(yr, ymin=low, ymax = up), fill = "orange", colour = "orange2") +
      geom_line(aes(yr,per), col = "orange3", size = 0.5 ) +
      {if (reverse)scale_x_reverse()} +
      scale_y_continuous(breaks =seq(0,100,10)) +
      aes(ymin=11) +
      xlab(ytitle) + ylab("%") +
      coord_flip() +
      theme_new +
      facet_grid(~forcats::fct_inorder(df$taxa),scales = "free", space = "free")

  }

  # stacked version
  if (type == "stacked") {

    # set default colors (if none were provided)
    cols = terrain.colors(nlevels(df$taxa))

    # and plot
    p <- ggplot(df, aes(x=yr, y=mean, fill = forcats::fct_inorder(taxa))) +
      geom_area() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      xlab(ytitle) + ylab("%") +
      theme_new
  }

  return(p)

}


#' plotMP
#'
#' Plots landcover reconstructed with MarcoPolo
#'
#' @param d the output of MarcoPolo application
#' @param ytitle  y-axis title, defaults to "Age cal. BP"
#'
#' @return A plot object
#'
#' @author Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
#' @author Almut Mrotzek <almut.mrotzek@gmx.de>
#' @references Mrotzek A., Couwenberg J., Theuerkauf M. and Joosten H. (2017) MARCO POLO –
#'   A new and simple tool for pollen-based stand-scale vegetation reconstruction. The Holocene 27: 321-330.
#'   DOI:10.1177/0959683616660171
#'
#' @examples ## Apply MarcoPolo on package data (from Mrotzek et al. 2017) and plot:
#' data("MP_example")
#'
#' ## run MarcoPolo
#' a <- MarcoPolo(localMP, regionalMP, paramsMP, assignmentMP)
#'
#' ## prepare plot
#' p <- plotMP(a, ytitle = "Site")
#'
#' # and write to pdf file
#' pdf("MarcoPolo plot.pdf", height = 12, width = 12)
#' p
#' dev.off()
#'
#' @export


plotMP <- function(d, ytitle = "Age cal. BP"){

  # plotting a stratigraphic diagramm with seperate colums for each taxon
  # inspired by code from Olivier Blarquez
  # from http://paleoecologie.umontreal.ca/684-2/

  # requires ggplot2 - is now in depends
  # if (!requireNamespace("ggplot2", quietly = TRUE)) {
  #   stop("Package ggplot2 needed for this function to work. Please install it.",
  #        call. = FALSE)
  # }

  # requires forcats
  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop("Package forcats needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #extracting the data
  d <- t(d)

  # get ages, if not numeric, just make sequence 1, 2, ...
  age <- rownames(d)
  #if (!is.numeric(age)) age <- seq(1:length(age))

  # now create data frames for median, 10% and 90% quantile
  # warum sortiert wer hier??
  df = data.frame(sample = rep(age,ncol(d)),
                  value = as.vector(as.matrix(d)),
                  taxa = as.factor(rep(colnames(d),each=nrow(d))))


  # define a theme without gridlines, etc. and elements that will decrease readability of the diagram
  theme_new <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # remove grids
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.x = element_text(size = 10, angle = 90, vjust = 0, hjust = 0.2), # Taxa names
                     strip.background = element_blank(),
                     strip.text.y = element_text(angle = 0),
                     legend.position = "right", panel.border = element_blank(),
                     axis.text.x = element_text(angle=0,hjust=0.5)) # Axis tick label angle


  # set default colors (if none were provided)
  cols = terrain.colors(nlevels(df$taxa))

  # and now the plot
  p <- ggplot() + geom_bar(aes(y = value, x = sample, fill = taxa), data = df,
                            stat="identity") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    xlab(ytitle) + ylab("%") +
    theme_new

  return(p)

}
