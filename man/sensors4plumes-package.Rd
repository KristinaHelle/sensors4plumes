\name{sensors4plumes-package}
\alias{sensors4plumes-package}
\alias{sensors4plumes}
\docType{package}
\title{
Test and optimise sampling designs based on plume simulations
}
\description{
The motivation of this package comes from radiological emergency preparedness, the functionality was needed for the planning of a network of sensors for the atmospheric gamma dose rate. These sensors should be placed in optimal locations to provide sufficient information in 'all' possible accident scenarios in the area of interest.
For this purpose, a model for the potential radioactive pollution is needed; as there is a wide variety of possible accidents, this model has to be stochastic. We know that pollution spreads by the physical laws of dispersion and transport, uncertainty comes from parameters like the wind field, the amount and kind of emitted pollutant and maybe the location of the source. It is difficult to find a parametric model for the distribution of possible pollutions, but we may describe it by a representative set of sample scenarios that can be simulated numerically, varying the uncertain parameters.
Each simulation represents a possible scenario for which we can check how well a proposed set of sensors can extract the required information: given the simulated values at the sensor locations we can mimic the process of information retrieval to decide if the result (trigger alarm, delineate evacuation area, etc.) fits the required result given the scenario. Quantifying this fitness and averaging it over the simulations, provides a global estimate about the fitness of a sensor set. This can guide algorithms to find optimal sampling designs.

The functionality of \code{sensors4plumes} is suited for the planning of monitoring sensor networks in all cases where sensors have to cope with a variety of possible scenarios that can be described by a set of simulations rather than by a parametric field.
The package provides functions to load such simulations automatically. To handle data that does not fit into memory, it resorts to the \link[raster]{raster-package}. It provides several basic functions to quantify the fitness of sensor sets that may be modified by users. These or user-defined fitness functions can arbitrarily be combined with various optimisation algorithms.
}
\details{
\tabular{ll}{
Package: \tab sensors4plumes\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2017-03-28\cr
License: \tab \cr
}
%The package can use simulation data both from GeoTiff and ascii files. If the data are too big to be loaded into memory, functionality of the \link{raster} package is used to handle it. A set of functions like \link{extractValues} help to load the data and turn it into an appropriate form for further computations.

%The values may refer to any spatial units: to points or grid cells, but also to polygons or lines.  These spatial objects are based on the \link{sp} package, extended by two new classes and generalised as \link{SpatialDataFrame}.
%The geometry may be given in the GeoTiff files or in separate objects. The spatial information and the simulated values are linked in \link{Simulations}.

%To evaluate sensor networks, the package provides a set of \link{costfunction}s. That basically summarise the simulated values at the sensor locations, maps generated from these data, or focus on their spatial relationship. Users may extend these functions or create their own.

%We provide a set of basic algorithms to \link{optimise} the sensor networks to have minimal cost. The resulting optima can be visualised or exported.

}
\author{
Kristina B. Helle

Maintainer: Kristina B. Helle <kristina.helle@uni-muenster.de>
}
%\references{
%~~ Literature or other references for background information ~~
%}
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
