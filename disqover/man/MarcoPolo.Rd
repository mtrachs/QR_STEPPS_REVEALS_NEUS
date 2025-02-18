% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MarcoPolo.R
\name{MarcoPolo}
\alias{MarcoPolo}
\title{MAnipulating pollen sums to ReCOnstruct POllen of Local Origin}
\usage{
MarcoPolo(local, regional, params, assignment)
}
\arguments{
\item{local}{local pollen data}

\item{regional}{regional pollen data}

\item{params}{table of R-values}

\item{assignment}{table naming the regional record for each local record}
}
\value{
The reconstructed local vegetation composition
}
\description{
MARCO POLO (MAnipulating pollen sums to ReCOnstruct POllen of Local Origin) is a quantitative approach
to stand-scale palynology, which translates pollen deposition from very small mires or ponds into
local vegetation composition. To that end, MARCO POLO removes those pollen types from the pollen sum
whose values are significantly higher than in a neighbouring large basin. The resulting regional pollen sum
is free of the disturbing factor of (extra-)local pollen. Based on this sum, comparison with the
pollen record from the large basin allows calculating sharp (extra-)local signals.
Treating the (extra-)local pollen portion with representation factors (R-values) then produces a
quantitative reconstruction of the stand-scale vegetation composition.
}
\examples{
  #### apply Marco Polo with the sample data from Mrotzek et al. (2017)
  data("MP_example", package = "disqover")
  a <- MarcoPolo(localMP, regionalMP, paramsMP, assignmentMP)

  #### data may be read from seperate files or from a single file with 4 sheets, as in the following example.
  #### for data structure see example data (e.g. MarcoPoloExampleData.xlsx)
  datafile <- system.file("extdata/marcopolo", "MarcoPoloExampleData.xlsx", package = "disqover")
  localMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="local")))
  regionalMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="regional")))
  paramsMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="params")))
  assignmentMP <- do.call(readxl::read_excel,args=c(list(datafile,sheet="assignment")))
  a <- MarcoPolo(localMP, regionalMP, paramsMP, assignmentMP)

  ### plot the results for surface samples
  ### prepare the plot
  p <- plotMP(a, ytitle = "Site")

  ### and write plot to a pdf file in the working directory
  pdf("MarcoPolo plot.pdf", height = 12, width = 12)
  p
  dev.off()


}
\references{
Mrotzek A., Couwenberg J., Theuerkauf M. and Joosten H. (2017) MARCO POLO –
  A new and simple tool for pollen-based stand-scale vegetation reconstruction. The Holocene 27: 321-330.
  DOI:10.1177/0959683616660171
}
\author{
Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>

Volkmar Liebscher <volkmar.liebscher@uni-greifswald.de>

Almut Mrotzek <almut.mrotzek@gmx.de>
}
