\name{interpol0reg}
\alias{interpol0reg}

\title{
Interpolate region-wise
}
\description{
\code{interpol0reg} is an interpolation function similar to and based on \code{\link[gstat]{idw0}} or \code{\link[gstat]{krige0}}. But the region of interest for the input as well as for the output may be split and interpolated with different interpolation methods. This can be useful to limit input to a subregion or to use different variograms. Altogether several interpolations may be applied and the results are combined.
It can be used as \code{fun_interpolation} in \code{\link{interpolate}}; therefore splitting is based on element names as explained below.
}
\usage{
interpol0reg(fun_interpolation = idw0z, y, data, newdata, 
  dataLoc, newdataLoc = 1:length(newdata), dataSplit = 1:length(newdata), newdataSplit = 1:length(newdata))
}
\arguments{
\item{fun_interpolation}{interpolation \code{function} (or list of those) with parameters \code{y}, \code{data}, \code{newdata} (form as for \code{\link[gstat]{krige0}}). All other parameters need default values that can be set by \code{replaceDefault} with \code{type = "interpolation_fun.interpolate"}. For \code{idw0z} see \code{\link{interpolate}}.}
\item{y}{\code{matrix} of input values, to be forwarded to \code{fun_interpolation}}
\item{data}{locations of input, best given as \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sp]{SpatialPixelsDataFrame}}, other possibilities described in \code{\link[gstat]{krige0}} may work. Number must agree with rows of \code{y}.}
\item{newdata}{locations of output, best given as \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sp]{SpatialPixelsDataFrame}}.}
\item{dataLoc}{\code{integer} names of input (\code{data}, \code{y}), length must fit input. Splitting by \code{dataSplit} refers to these names. - This can be used to allow hierarchical selection (as inside of \code{\link{interpolate}}) and then applying the splitting to the selection.}
\item{newdataLoc}{\code{integer} names of output (\code{newdata}), length must fit. Splitting by \code{newdataSplit} refers to these names.}
\item{dataSplit}{list of \code{integer} names of input to be used in the different interpolations.}
\item{newdataSplit}{list of \code{integer} names of output to be generated in the different interpolations.}
}

\details{
\code{newdataSplit} splits the set of locations where interpolation is to be computed. For each list entry one interpolation method is applied. The function and the indices of the input for each interpolation are taken from the respective list entries in \code{fun_interpolation} and \code{dataSplit}. If these lists do not have the same length, an error is thrown. Instead of lists you may use one function or one list of indices instead which then is used in all interpolations. The results of all interpolations are combined to return values for all locations in \code{newdata}. If output names from \code{newdataLoc} do not occur in any entry of \code{newdataSplit} \code{NA} is returned there, if they occur in several, the last interpolation determines the returned value (with warning). To use it inside of \code{\link{interpolate}}, all parameters except \code{y, data, newdata, dataLoc, newdataLoc} must be set beforehand by \code{\link{replaceDefault}} using \code{type = "fun_interpolationSplit.interpolate"}
}

\value{
A matrix of interpolated values at the locations of \code{newdata} (rows); each column generated from the respective column in \code{y}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{ kristina.helle@uni-muenster.de  }

\examples{
# get data
data(SPointsDF)

# define names and splitting
# input
d = c(1,2,3,8,9)
dS = list(dS1 = 7:9, dS2 = 1:6)
# output
n = 1:8
nS = list(nS1 = 1:3, nS2 = c(2,5:9)) # overlapping and not covering all elements 
# functions (kriging with different variograms)
fun1 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(10, "Exp", 50000, 5), ... = NA), type = "fun_interpolation.interpolate")[[1]]
fun2 = replaceDefault(krige0, newDefaults = list(formula = z ~ 1, beta = NA, model = vgm(100, "Sph", 500, 50), ... = NA), type = "fun_interpolation.interpolate")[[1]]

# use different input and function in the 2 interpolations
int_1 = interpol0reg(fun_interpolation = list(fun1,fun2),
                         dataSplit = dS,
                         newdataSplit = nS,
                         y = as.matrix(SPointsDF@data)[d,, drop = FALSE], 
                         data = SPointsDF[d,], 
                         newdata = SPointsDF[n,], 
                         dataLoc = d, 
                         newdataLoc = n) 

# use the same input and different functions                         
int_2 = interpol0reg(fun_interpolation = list(fun1,fun2),
                         dataSplit = d,
                         newdataSplit = nS,
                         y = as.matrix(SPointsDF@data)[d,, drop = FALSE], 
                         data = SPointsDF[d,], 
                         newdata = SPointsDF[n,], 
                         dataLoc = d, 
                         newdataLoc = n) 
}                         