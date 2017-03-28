\name{similaritySD}
\alias{similaritySD}

\title{
Determine (spatial) similarity between sampling designs
}
\description{
Takes sampling designs and a reference sampling design and computes similarity between them by various algorithms: number of coinciding sensors, sensors within a k-neighbourhood (based on \code{\link[FNN]{get.knnx}}), Earth Mover's Distance (based on \code{\link[emdist]{emd2d}}).
}
\usage{
similaritySD(simulations, SD, referenceSD, type = "equal", k = 9, ...)
}

\arguments{
  \item{simulations}{\code{\link{Simulations}} object, needed to define spatial properties; not needed if \code{type = "equal"}}
  \item{SD}{vector, matrix (sampling designs in rows) or list of sampling designs; indices of the \code{locations} of \code{simulations}}
  \item{referenceSD}{sampling design (indices of the \code{locations} of \code{simulations}); all \code{SD} are compared to this sampling design}
  \item{type}{character, indicating comparison method, one of \code{"equal"}, \code{"Kneighbours"}, \code{"EarthMoversDistance"}}
  \item{k}{neighbourhood size, to be used for \code{type = "Kneighbours"}, forwarded to \code{\link[FNN]{get.knnx}}}
  \item{...}{parameters to be used for \code{type = "EarthMoversDistance"}, forwarded to \code{\link[emdist]{emd2d}}}
}
\details{
Originally Earth Mover's Distance compares matrices of two distributions. We can use it for cases where all locations lie on a regular grid. Then a sampling design is represented by a matrix of all grid cells: each sensor is represented by 1, all other cells are 0 (for locations on a \code{SpatialPolygrid}, a sensor may belong to several cells, then the weight 1 is distributed evenly to all these cells).
}
\value{
A vector of length referring to the number of sampling designs in \code{SD} with the similarity between these sampling designs and \code{referenceSD} respectively.
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
## ---------------- prepare data --------------------------
# simulations with SpatialPolygridDataFrame
data(radioactivePlumes)

# simulations with SpatialPixelsDataFrame
data(meuse.grid)
coordinates(meuse.grid) = ~ x + y
gridded(meuse.grid) = TRUE
meuseGSim = subset(radioactivePlumes, 
                  locations = 1:length(meuse.grid), plumes = 1:10)
meuseGSim@locations = meuse.grid

# simulations with rectangular SpatialPixelsDataFrame
data(SPixelsDF)
rectSim = subset(radioactivePlumes, 
                 locations = 1:length(SPixelsDF), plumes = 1:10)
rectSim@locations = SPixelsDF # not nice, no checking if sizes fit

# simulations with SpatialPointsDataFrame 
data(meuse)
coordinates(meuse) = ~ x + y
meuseSim = subset(radioactivePlumes, 
                  locations = length(meuse), plumes = 1:10)
meuseSim@locations = meuse # not nice, no checking if sizes fit

# for meuseSim
sds0 = c(10,  25,  42,  84,  90,  92,  94,  97, 120, 153)#sample.int(155, 10)
sds1 = matrix(c(97,   79,   40,   68,  131,
                57,   18,   38,  118,   14,
                23,   71,   22,   94,   27,
                125,  108,    4,   80,  129,
                130,   96,  101,  137,   19,
                77,  138,   32,   95,   88,
                140,   73,   43,  153,   12,
                8,  141,   92,   35,  102), 
              byrow = TRUE, nrow = 8)#matrix(sample.int(155, 40), nrow = 8)
sds2 = list(
  c(28, 59, 64, 78, 81),
  c(5,  13,  21,  31,  45,  91,  92, 122, 130, 141),
  c(1,   2,  19,  36,  50,  51,  58,  59,  90, 103, 105, 107, 112, 123, 132),
  c(34,  48,  49,  50,  71,  76,  77,  86,  92,  97, 100, 103, 104, 106, 108, 
  110, 113, 127, 134, 142)
)

# ------------- "equal" -----------------------------
  sim_1_1 = similaritySD(# no simulations needed
    SD = sds1, referenceSD = sds0, type = "equal")  

#------------- "Kneighbours" -------------------
  sim_2_1_1 = similaritySD(simulations = meuseGSim,
                              SD = sds2, referenceSD = sds0, 
                              type = "Kneighbours", k = 4)  
  sim_2_1_2 = similaritySD(simulations = meuseGSim,
                         SD = sds2, referenceSD = sds0, 
                         type = "Kneighbours", k = 9) 

# plot points for visual comparison
plot(coordinates(meuseSim@locations), pch = ".")
   points(coordinates(meuseSim@locations)[sds0,], cex = 0.5)
   for (i in seq(along = sds2)){
     points(coordinates(meuseSim@locations)[sds2[[i]],], 
            col = i + 1)
   }
#---------- "EarthMoversDistance"----------
  # SpatialPixelsDataFrame
  sim_3_1_3 = similaritySD(simulations = rectSim,
                                SD = c(2,9), referenceSD = c(3,8), 
                                type = "EarthMoversDistance")

  sim_3_1_4 = similaritySD(simulations = rectSim,
                                SD = c(3,8), referenceSD = c(4,5), 
                                type = "EarthMoversDistance")
  sim_3_1_5 = similaritySD(simulations = rectSim,
                                SD = c(3,8), referenceSD = c(4,5), 
                                type = "EarthMoversDistance", dist = "manhattan")  

  ## SpatialPolygridDataFrame
  sim_3_1_5 = similaritySD(simulations = radioactivePlumes,
                                SD = c(220, 503), referenceSD = c(224, 544), 
                                type = "EarthMoversDistance")  

}
