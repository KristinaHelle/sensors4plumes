\name{loadSimulations}

\alias{loadSimulations}

\title{
Load values from raster or text files into Simulations objects
}

\description{
Plume simulations for one scenario are a set of maps of the same area, each showing one plume. 
This function generates a \code{Simulations} object from all data of a scenario which can then be used to compute cost and to optimise sampling designs.
Data may extend the size of the memory; it can be provided either as raster maps or as values in text files (in this case spatial reference can be added later). The resulting \code{values} are saved to disk as a \code{\link[raster]{raster}} object.
}

\usage{
loadSimulations(basicPath = ".", readBy = NA, 
  multilayer = "kinds", region, bBox, 
  nameSave = NA, overwrite = FALSE, ...)
%loadSimulations(basicPath = ".", 
%  readBy = c(NA, "scan", "raster"), 
%  multilayer = c("kinds", "plumes"), 
%  region, bBox, 
%  savePath = ".", overwrite = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{basicPath}{path to directory with files; it must not contain any but the files with simulations to be used
  }
  \item{multilayer}{way to interpret files with several layers/columns; \code{"kinds"} (default): layers/columns are interpreted as different kinds of values (number of layers/columns must be the same in all files) or \code{"plumes"}: layers/columns are interpreted as values of different plumes, all being of the same kind
  }
  \item{readBy}{method to read the files; for some suffixes the reading method is determined automatically (all files must have the same suffix): \code{"scan"} for \code{".txt"} or \code{".csv"}, and \code{"raster"} for \code{".tif"} or \code{".grd"}; else it has to be specified here, (default is \code{NA})  
  }
  \item{region}{logical indicating in which locations to keep the values (for \code{raster} files, order is row-wise), length must equal number of locations (if combined with \code{bBox} length must fit locations inside \code{bBox})
  }
  \item{bBox}{bounding box of locations to keep, given as vector \code{c(xmin, xmax, ymin, ymax)}; works only for raster data, for text files it is ignored
  }
  \item{nameSave}{name to save resulting raster files that are created automatically and then contain the final data, files are created by appending names to it, starting with \code{"_"}
  }
  \item{overwrite}{logical if files at \code{nameSave} may be overwritten
  }
  \item{...}{parameters to be forwarded to read.table and scan like \code{skip}, \code{sep}, \code{dec}}
}

\details{
Plume simulations have to be a set of maps all defined by values at the same set of locations; there may be several kinds of maps for each plume. 
Such data can be loaded from files automatically if these have the following properties:\cr
Files can either all be loaded by \code{\link[raster]{brick}} (.grd/.gri, .tif) or by \code{\link{scan}} (.txt, .csv). If they have one of the listed suffixes, the method is derived automatically; for other suffixes it must be given by the user in \code{readBy}.\cr
Each file must contain the values of a plume in all locations. \cr
Files may contain several kinds of values; in this case each layer (of .tif or .grd files) or each column (in .txt or .csv files) belongs to different kind of map, the structure (number and meaning of columns/layers and rows) must be common to all files. \cr
If there is only one kind of map, files may contain columns/layers with the values of several plumes.\cr
File names must follow some rules: 
They must not contain any \code{"."} except the one to start the suffix. 
They must not contain \code{"_"}, except for the case described below; they are interpreted as the names of plumes. 
If there are several kinds of maps, either all kinds for a plume are in the same file or there is a separate file for each kind and each plume. In the latter case file names must be of the form \code{plumeName_kindName.suffix}. For each plume there must be a file for each kind. \cr
}

\value{
The output of this function depends on the file type: raster files result in a \code{Simulations} object, with spatial properties kept in the \code{locations} (projection etc. of GeoTiff); text files do not contain spatial information, therefore a \code{list} of \code{plumes} and \code{values} is returned that can be used to generate a \code{Simulations} object together with appropriate \code{locations}.

The function may generate (and overwrite or delete) files at \code{nameSave}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
# examples for possible input is given in package sensors4plumesData in the 
# subfolders of inst/extdata/fileFormats

# from raster files
Simulations1 = loadSimulations(
  basicPath = paste0(path.package("sensors4plumesData"),
                     "/extdata/fileFormats/raster"),
  nameSave = "t1_",
  overwrite = TRUE
)

# from text files
Simulations2 = loadSimulations(
  basicPath = paste0(path.package("sensors4plumesData"),
                     "/extdata/fileFormats/text_multicolumn"),
  readBy = "scan",  
  nameSave = "t2_",
  overwrite = TRUE,
  skip = 1
)

# from multilayer raster files, selecting locations in a bounding box
region = as.logical(replicate(20, sample.int(2,1)) - 1)
bBox = c(0.2, 0.7, 0.4, 0.8)

Simulations3 = loadSimulations(
  basicPath = paste0(path.package("sensors4plumesData"),
                     "/extdata/fileFormats/raster_multilayer"),
  nameSave = "t3_",
  region = region,
  bBox = bBox,
  overwrite = TRUE
)
}
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
