% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotREVEALS}
\alias{plotREVEALS}
\title{plotREVEALS}
\usage{
plotREVEALS(d, type = "single", ytitle = "Age cal. BP", reverse = TRUE,
  error = "quantiles")
}
\arguments{
\item{d}{the output of REVEALSinR}

\item{type}{plot type: 'single' for a stratigraphic diagram with seperate
columns. Error plotted as 10\% and 90\% qantiles; 'stacked' for stacked area plot.
Colours by default taken from the terrain.colors palette}

\item{ytitle}{y-axis title, defaults to "Age cal. BP"}
}
\value{
A plot object
}
\description{
Plots landcover reconstructed with REVEALSinR.
}
\examples{
## Apply REVEALSinR and plotting on a core section from Lake Tiefer See:

# load 'Tiefer See' data set and apply REVEALSinR
data("Tiefer_See")
cover <- REVEALSinR(pollenTS, paramsTS, tBasin = "lake", dBasin = 600)

# produce stratigraphic plot (like a pollen diagram)
# by default, the y-axis is reversed to suit BP dates
p <- plotREVEALS(d = output, type = "single", ytitle = "Age CE", reverse = FALSE)

# and write to pdf file
pdf("reveals_single.pdf", height = 8, width = 12)
p
dev.off()

# now produce stacked area plot
p <- plotREVEALS(d = output, type = "stacked", ytitle = "Age CE")

# and write to pdf file
pdf("reveals_stacked.pdf", height = 12, width = 12)
p
dev.off()

}
\author{
Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>

Triin Reitalu <triinreitalu@gmail.com>
}
