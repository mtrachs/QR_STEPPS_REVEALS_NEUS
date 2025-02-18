% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edar.R
\name{edaimportxl}
\alias{edaimportxl}
\title{Preparing data for the Extended Downscaling Approach with edar()}
\usage{
edaimportxl(pollenfile, envfile, paramsfile)
}
\arguments{
\item{pollenfile}{path and name of the pollen file. This file includes one
sheet with raw pollen counts for each taxon and site (usuall lake). Each
row represents one pollen taxon, each column one site. The first row
includes site names, the first column taxon names. See
'eda_tc17c_pollen.xlsx' as example.}

\item{envfile}{path and name of the environmental data file. This file
includes spatial information, i.e. the cover of each environmental
unit/class (e.g. each soil type)  within a specific distance from each
site. The area is given as the total area inside the respective full circle
around each site (not in donuts). The file includes one sheet for each
environmental unit. In each sheet, each row represents one radius (in
increasing order), each column one site. The first row includes site names,
the first column radii (as integer or real number). Sheet names should
represent the environmental unit/class. See 'eda_tc17c_environment.xlsx' as
example.}

\item{paramsfile}{path and name of the parameter file. This file includes
fallspeed of pollen (in m/s), pollen productivity estimates and standard
error of the PPEs. PPEs are region specific corrections factors, so values
from a suitable calibration study are required. PPE calibration is highly
sensitive to dispersal model selection. Therefore, the same underlying
dispersal model has to be applied in PPE calibration and EDA application.
The parameters have to be provided in three columns: 'fallspeed', 'ppe'
and 'ppe.error'. Names of the pollen taxa are represented in row names. See
'eda_tc17c_params.xlsx' as an example. The parameter file may include more
pollen taxa than the pollen data, and in a different order.}
}
\value{
...edav: a list with input data for edar()
}
\description{
Applying the Extended Downscaling Approach (EDA) with the edar() function
requires a number of data: pollen data, environmental data and taxon
parameters. All the data has to be provided in a single list. edaimportxl
builds this list from three data files. Two of them, the file for pollen data
and the file for taxon parameters shall only include one sheet. The file for
environmental data includes a number of sheets, one for each class/unit (e.g.
soil type). For the structure of each file, see the examples from Theuerkauf
& Couwenberg (2017), which are provided as external data.
}
\examples{
# edaimportxl() prepares data for edar()
# it requires three files: pollen data, environmental data and taxon parameters
# the following example uses external data included in the package

# first set path and filenames for the data files
pollenfilename <- system.file("extdata/eda", "eda_tc17c_pollen.xlsx", package = "disqover")
envfilename <- system.file("extdata/eda", "eda_tc17c_environment.xlsx", package = "disqover")
paramsfilename <- system.file("extdata/eda", "eda_tc17c_params.xlsx", package = "disqover")

# now use edaimportxl()
eda_tc17c <- edaimportxl(pollenfilename, envfilename, paramsfilename)

}
\references{
Theuerkauf M. and Couwenberg J. (2017) The extended downscaling
  approach: A new R-tool for pollen-based reconstruction of vegetation
  patterns. The Holocene 27(8): 1252–1258.
  https://doi.org/10.1177/0959683616683256
}
\author{
Martin Theuerkauf <martin.theuerkauf@uni-greifswald.de>
}
