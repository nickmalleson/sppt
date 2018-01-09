#' @import sppt
NULL


#
# NOTE:
# The expanding cell code has originally come from the project: spatialtest/r/expanding_cell.Rmd
# XX LINK TO GITHUB REPO
# That was under development, so I have copied all of the code over and adapted it here.
#

#' An extension to the sppt package that runs the test over a number of grids with varying
#' resolutions.
#' TODO - integrate this script into the sppt package properly.
#'
#' @param resolution.cutoff The maximum number of cells to use. I.e. in the results there will not be
#' a grid with more than 'resolution.cutoff' cells.
#'
#' @return A list of lists. The first list dimension is for each resolution, and the second dimension
#' has the individual grids that make up the results at that resolution.

sppt.exp.cell <- function(base_points.sp, test_points.sp, uoa.sp, resolution.cutoff=1000,
                          outputlist=FALSE, nsamples=200, percpoints=85, conf_level=95){

  # Create the regular grid at the finest resolution (and make it slightly larger than necessary)

  XXXX

  # Create four additional grids by shifting N, S, E, W

  XXXX

  # Run the test on those four grids

  XXXX

  # Iteratively aggregate the grids until the grid only has one cell in it

  XXXX

  #  Return the results list.
}



# Stuff to test the algorithm


library(GISTools)
library(rgdal)     # For reading shapefiles
library(raster)    # For creating regular grids
library(classInt) # Jenks natural breaks
library(hydroGOF)   # Has an rmse() function

setwd("/Users/nick/research_not_syncd/git_projects/sppt")



points1 <- readOGR(dsn = "./data-raw/basic_test/", layer = "points1")
points2 <- readOGR(dsn = "./data-raw/basic_test/", layer = "points2")

plot(points1, col="blue")
points(points2, col="red")

# Create the regular grid at the largest resolution (and make it slightly larger than necessary)

bb <- bbox(points1 + points1) # A bounding box around all points


XXXX

# Create four additional grids by shifting N, S, E, W

XXXX

# Run the test on those four grids

XXXX

# Iteratively make the grids smaller XXXX the grids until the grid only has one cell in it

XXXX

#  Return the results list.



