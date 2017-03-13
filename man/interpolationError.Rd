\name{interpolationError}
\alias{interpolationError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate many maps at once, compare them to the original and determine a global error
}
\description{
This function calls \code{interpolate}, then it compares the result plume-and-location-wise
to the original and summarises the resulting error values.
}
\usage{
interpolationError(simulations, locations, values, 
  fun_interpolation, fun_error, fun_Rpl, fun_Rpl_cellStats = "mean", fun_l, 
  tmpfile = "tmp_interpolationError", overwrite = FALSE, chunksize = 1e+7)
}
\arguments{
  \item{simulations}{\code{Simulations}}
  \item{locations}{indices of \code{locations} to be used as input; multiple and invalid values are ignored}
  \item{values}{ layer of the \code{values} of \code{simulations} to be used; interpolation can only be applied to one layer: if \code{values} is a vector, only the first entry is used }
  \item{fun_interpolation}{ interpolation \code{function}, must have the parameters \code{y}, \code{data}, \code{newdata} (form as for \code{\link[gstat]{krige0}}). All other parameters need default values, e.g. the \code{model} of \code{krige0} that can be set by \code{replaceDefault} with \code{type = "interpolation_fun.interpolate"}. }% output??
  \item{fun_error}{ function to compare original and interpolated map location-and-plume-wise; must have a parameter \code{x}, then \code{x[1]} is the original and \code{x[2]} the interpolated value; it is forwarded to \code{\link{simulationsApply}} as \code{fun_pl}, therefore it has to fulfil all requirements for such functions; \code{\link{interpolationErrorFunctions}} provides some common examples}
  \item{fun_Rpl}{ function to summarise the location-and-plume-wise errors, is forwarded to \code{simulationsApply} as \code{fun_Rpl}; if input does not fit into memory, it cannot be applied (causes warning) }
  \item{fun_Rpl_cellStats}{ alternative function to summarise the location-and-plume-wise errors, is forwarded to \code{simulationsApply} as \code{fun_Rpl_cellStats}}
  \item{fun_l}{ function to compre original and interpolated location-wise, i.e. generate one global map that takes into account all plumes; must have parameter \code{x}, then \code{x[,1]} refers to the original and \code{x[,2]} to the interpolated values; it is forwarded to \code{\link{simulationsApply}} as \code{fun_pl}, therefore it has to fulfil all requirements for such functions; \code{\link{interpolationErrorFunctions}} provides some common examples }
  \item{tmpfile}{ filename for the raster file in case the result does not fit into memory; if \code{FALSE} the function stops with a warning and does not create a file }
\item{overwrite}{ boolean, if the file at \code{tmpfile} may be overwritten } 
\item{chunksize}{ maximal number of cells to be processed at once -- forwarded to \code{\link[raster]{blockSize}} inside }
}

\value{
List of values and rasters (of same dimension as the \code{values} of the \code{simulations}):

  \code{"cost"}: result of \code{fun_Rpl} if available (if not, warning), else result of \code{fun_Rpl_cellStats} (to guarantee that there is always a value)
  
  \code{"cost_cellStats"}: result of \code{fun_Rpl_cellStats} (if this is not in \code{"cost"}) 
  
  \code{"error_locationsplumes"}: raster, result of \code{fun_error} 
  
  \code{"interpolated"}: result of the interpolation with \code{fun_interpolation}
  
  \code{"costLocations"}: result of \code{fun_l}

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{ kristina.helle@uni-muenster.de  }

\examples{
  data(radioactivePlumes_area)
  ## preparation
  idw0z = replaceDefault(idw0, newDefaults = list(
    formula = z ~ 1))[[1]]
  sampleLocations250 = sample.int(nLocations(radioactivePlumes_area), 250)
  fun_Rpl_mean = function(x, nout = 1){ 
    mean(x[,1], na.rm = TRUE)
  }
  ## interpolation error (computation in memory)
  interpolationError_delineation <- interpolationError(
      simulations = radioactivePlumes_area,
      locations = sampleLocations250,
      values = 2,
      fun_interpolation = idw0z,
      fun_error = delineationError,
      fun_Rpl = fun_Rpl_mean,
      fun_Rpl_cellStats = "mean",
      fun_l = delineationErrorMap
    )    
  ## plot error map
  interpolationErrorMaps = radioactivePlumes_area
  interpolationErrorMaps@values = stack(radioactivePlumes_area@values[[2]],
                                        interpolationError_delineation[["interpolated"]],
                                        interpolationError_delineation[["error_locationsplumes"]])
  interpolationErrorMapsSDF = extractSpatialDataFrame(interpolationErrorMaps, plumes = 1:5) 
  interpolationErrorMapsSDF@data$costMap = interpolationError_delineation[["costMap"]]
  spplot(interpolationErrorMapsSDF, zcol = c(1:10, 21:25))
  spplot(interpolationErrorMapsSDF, zcol = "costLocations")
}

