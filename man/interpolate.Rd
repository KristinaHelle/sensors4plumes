\name{interpolate}
\alias{interpolate}
\alias{fitMedianVariogram}
\alias{idw0z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate many maps at once -- even if they do not fit into memory
}
\description{
\code{interpolate} is a wrapper for \code{\link[gstat]{idw0}} and \code{\link[gstat]{krige0}} to interpolate several maps at once if locations of input values and desired output agree, even if output and maybe even input does not fit into memory. Works by writing output (and maybe input and intermediate results) to raster files.

\code{fitMedianVariogram} is a wrapper of \code{\link[automap]{autofitVariogram}}: fits variograms to a sample of plumes of a \code{simulations} object and generates a variogram with the median parameters.
}
\usage{
interpolate(simulations, locations, kinds = 1, fun_interpolation, 
  tmpfile = "tmp_interpolate", overwrite = FALSE, chunksize = 1e+7)
fitMedianVariogram(simulations, plumes, locations, kinds = 1)
idw0z(formula = z ~ 1, data, newdata, y, idp = 2)
}
\arguments{
  \item{simulations}{\code{Simulations}}
  \item{locations}{indices of \code{locations} to be used as input; multiple and invalid values are ignored}
  \item{plumes}{indices of \code{plumes} to fit variograms to}
  \item{kinds}{ layer of the \code{values} of \code{simulations} to be used; interpolation can only be applied to one layer: if \code{values} is a vector, only the first entry is used in both functions. }
  \item{fun_interpolation}{ interpolation \code{function}, must have the parameters \code{y}, \code{data}, \code{newdata} (form as for \code{\link[gstat]{krige0}}). All other parameters need default values, e.g. the \code{model} of \code{krige0} that can be set by \code{replaceDefault} with \cr \code{type = "interpolation_fun.interpolate"}. 
  
  %It may also have the two additional parameters \code{dataLoc} and \code{newdataLoc} (like \code{\link{interpol0reg}}) which can by modified in \code{\link{replaceDefault}} with \code{type = "fun_interpolationSplit.interpolate"}.
  }
  \item{tmpfile}{ filename for the raster file in case the result does not fit into memory; if \code{FALSE} the function stops with a warning and does not create a file }
\item{overwrite}{ boolean, if the file at \code{tmpfile} may be overwritten } 
\item{chunksize}{ maximal number of cells to be processed at once -- forwarded to \code{\link[raster]{blockSize}} inside }
\item{formula}{ formula that defines the dependent variable as linear model of independent variables, forwarded to \code{\link[gstat]{idw0}} -- see also there, default is no independent variables}
\item{data}{data frame with dependent variable and coordinates, forwarded to \code{\link[gstat]{idw0}} -- see also there}
\item{newdata}{data frame or Spatial object with prediction locations, forwarded to \code{\link[gstat]{idw0}} -- see also there}
\item{y}{matrix, forwarded to \code{\link[gstat]{idw0}} -- see also there}
\item{idp}{, forwarded to \code{\link[gstat]{idw0}} -- see also there}
}

\value{
\code{interpolate} returns a \code{\link[raster]{RasterLayer-class}} with the interpolations; it has the same size as the \code{values} of the input \code{simulations} and belongs to the same locations and plumes. Also projection is the same, so they can be combined by \code{\link[raster]{stack}}. If it does not fit into memory it is saved at \code{tmpfile} with the extension \code{"_interpolated.grd"}.
The function may produce intermediate files at \code{tmpfile} (with name extensions) that are deleted in the end.

\code{fitMedianVariogram} returns a variogram model (\code{\link[gstat]{vgm}}).

\code{idw0z} is \code{\link[gstat]{idw0}} with one extra default parameter: \code{formula = z ~ 1}. This way it can directly be used as \code{fun_interpolation} in \code{interpolate}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
\dontrun{
## takes some time
# get data
data(radioactivePlumes)

# generate median variogram from plumes
\dontrun{
## takes some seconds
medianVariogram = fitMedianVariogram(simulations = radioactivePlumes, 
                              plumes = 1:nPlumes(radioactivePlumes),
                              kinds = 1)
}
## the result is in:
data(medianVariogram)

# prepare interpolation function 
krige0var = replaceDefault(krige0, newDefaults = list(
  formula = z ~ 1, model = medianVariogram, beta = NA, ... = NA),
  type = "fun_interpolation.interpolate")[[1]]

# sample locations: proposed sensors
sampleLocations = sample.int(nLocations(radioactivePlumes), 50)

# interpolate
interpolated = interpolate(
  simulations = radioactivePlumes,
  kinds = 1,
  locations = sampleLocations,
  fun_interpolation = krige0var)
  
# combine plot original and interpolated
originalAndInterpolated = radioactivePlumes
originalAndInterpolated@values = stack(
  originalAndInterpolated@values[[1]], interpolated)

OriginalAndInterpolated = extractSpatialDataFrame(
  originalAndInterpolated, plumes = 1:4)

samplePoints = 
  as(OriginalAndInterpolated, "SpatialPointsDataFrame")[sampleLocations,]
spplotLog(OriginalAndInterpolated,
  sp.layout = list("sp.points", 
    samplePoints, col = 3))
spplot(OriginalAndInterpolated,
  sp.layout = list("sp.points", 
    samplePoints, col = 3))

}
}
