WORKING_DIR <- '~/mapping/projects/sppt/'
setwd(WORKING_DIR)

library(GISTools)
#library(rgeos)    # For things like gIntersects
library(rgdal)     # For reading shapefiles
library(raster)    # For creating regular grids
#library(plyr)     # For counting and aggregating
#library(tmap)     # For thematic maps
library(classInt) # Jenks natural breaks
#library(png)      # For loading pngs after they have been written
#library(grid)     # so that they can be embedded in knitted documents
#library(spdep)    # For doing the spatial regression, contiguity matrices, etc.
#library(GWmodel)  # For geographically weighted regression (GWR)
#library(MASS)     # For stepwise regression (stepAIC())
#library(pander)    # For printing tables nicely
#library(MVN)      # For testing for multivariate normality (MVN)
#library(RColorBrewer) # For making nice colour themes
library(hydroGOF)   # Has an rmse() function
library(xtable)   # For making latex/html tables
library(parallel) # For ruonning things in parallel (e.g. mclapply())
no_cores <- detectCores() # Detect the number of cores that are available and use half (often CPUs simulate 2 threads per core)
Sys.setenv(MC_CORES=no_cores) # Run on n cores (I'm not sure which of these
options("mc.cores"=no_cores) # is correct).
library(Deriv)   # For calculating deriviatives

library(sppt) # The spatial point pattern test library (see https://github.com/wsteenbeek/sppt for install instructions)

ZIPFILENAME <- "./data-raw/vancouver_public/crime_shp_all_years.zip"


# Unzip the files into the working directory
zipfile <- unzip(ZIPFILENAME)

# Read the shapefile
all.crime <- readOGR(dsn="./crime_shp_all_years", layer = "crime_shp_all_years")

# Delete the extracted file
unlink(zipfile)
unlink("crime_shp_all_years", recursive = TRUE) # a left over directory

# Assign short codes to crime types (the below throws an error about 'values must be length 1'
# if there are any crime types that aren't matched to short equivalents)
all.crime$TYPE2 <- as.factor(unlist(
  mclapply(X = all.crime$TYPE, FUN = function(t) {
    t2 <- switch(as.character(t),
          "Break and Enter Commercial" = "BNEC", 
          "Break and Enter Residential/Other" = "BNER",
          "Mischief" = "MISCHIEF",
          "Other Theft" = "OTHERTHEFT",
          "Theft from Vehicle" = "TFV",
          "Theft of Bicycle" = "TOB",
          "Theft of Vehicle" = "TOV",
          "Vehicle Collision or Pedestrian Struck (with Fatality)" = "COLFAT", 
          "Vehicle Collision or Pedestrian Struck (with Injury)" = "COLINJ")
})))
  
# Drop 2018 (not sufficient data yet)
all.crime <- all.crime[which(all.crime$YEAR<2018),]

rm(ZIPFILENAME)


# A directory for the crime data (in case it doesn't exist)
dir.create(file.path("./ROUT", "crime"), showWarnings = FALSE)

# Shorter versions of the crime types
crime.types2 <- unique(all.crime$TYPE2)

for (type in crime.types2) {
  assign(as.character(type), all.crime[all.crime$TYPE2==type,])
  # Write them as well as they can be useful later
  writeOGR(all.crime[all.crime$TYPE2==type,], dsn = "./ROUT/crime", layer = type, 
           driver = "ESRI Shapefile", overwrite_layer = TRUE)
}






msea <- function(points1, points2, N=20, n.shifts = 10, ignore.zeros=FALSE, step=1, 
                 return.sobject=TRUE, return.grids=TRUE) {
  
  # Check that the projections are the same
  if ( proj4string(points1) != proj4string(points2) ) {
    warning("The points1 and points2 projections are different, this will probably lead to catastrophic results!")
    stop()
  }
    
  bb <- bbox(points1 + points1) # A bounding box around all points
  
  # Store all grids (data frames) in a big long list
  results <- list()
  # Remember some other things that are useful later
  cell.areas <- c() # The area of the cells
  num.cells <- c()  # The number of cells in each iteration

  # Remember the global errors associated with each grid
  rss <- c() # Residual sum of squares
  r.squared <- c()
  rmse <- c()
  globalS <- c()
  globalS.robust <- c()
  
  # Keep a link to the object that is returned from the call to sppt. Useful for debugging mostly.
  s.object <- c()
  
  # Remember the iteration and shift numbers (these are the i and j in the nested loops)
  iteration <- c()
  shift <- c()
  
  # Create the grids - adapted from Brunsdon & Comber (2015, p150)
  # Note: will actually start at i=2 which gives 2*2=4 cells (doesn't make sense to calculate error for 1 data point)
  # but it's easier to start from i=1 and delete that result afterwards
  counter <- 1 # for counting the total number of grids created (required for indexing)
  for (i in seq(from=1,to=N,by=step)) {
    # Cell size is the total width divided by the number of cells to draw so far (i)
    cell.width <-  (bb[1,2] - bb[1,1]) / i
    cell.height <- (bb[2,2] - bb[2,1]) / i
    
    # Make the bounding box slightly larger than necessary (by half a cell in each direction), 
    # so when the grid is shifted there wont be any points outside
    # It needs to be big enough so that it can have one extra ring of cells around it
    bb.larger <- bb # The new bounding box
    bb.larger["coords.x1","min"] <- bb["coords.x1","min"] - ( cell.width / 2  ) # Min x gets smaller
    bb.larger["coords.x2","min"] <- bb["coords.x2","min"] - ( cell.height / 2 ) # Min y gets smaller
    bb.larger["coords.x1","max"] <- bb["coords.x1","max"] + ( cell.width / 2  ) # Max x gets larger
    bb.larger["coords.x2","max"] <- bb["coords.x2","max"] + ( cell.height / 2 ) # Max y gets larger
    
    # For each resolution, repeat a few times by slightly shifting the grid by a random amount in a random direction
    for (j in 1:n.shifts) {
      # Remember the cell area (useful later) (needs to be repeated for each shift)
      cell.areas <- c(cell.areas, (cell.width * cell.height) )

      # Chose random N-S and E-W directions to shift the grid in (using a random uniform distribution)
      shift.x <- runif(n=1, min=-cell.width /2, max=cell.width /2 )
      shift.y <- runif(n=1, min=-cell.height/2, max=cell.height/2)
      
      # Calculate the centre of the lower-left cell (the one with the smallest coordinates),
      # taking into account the shift
      centre.x <- bb.larger[1,1] + ( cell.width  / 2 ) + shift.x
      centre.y <- bb.larger[2,1] + ( cell.height / 2 ) + shift.y
      
      # Create a grid  
      grd <- GridTopology(
        cellcentre.offset = c(centre.x, centre.y), # No offset, the grid will just cover all the points
        cellsize = c(cell.width, cell.height),
        cells.dim = c(i+1,i+1)
      )
      
      number.of.cells <- (i+1) * (i+1) # Add an extra row and column to account for shifting 
      num.cells <- c(num.cells, number.of.cells) # Remember the number of cells in this iteration
      
      # Convert the grid into a SpatialPolygonsDataFrame
      spdf <- SpatialPolygonsDataFrame(
        as.SpatialPolygons.GridTopology(grd),
        data = data.frame(c(1:number.of.cells)),
        match.ID = FALSE
      )
      proj4string(spdf) <- proj4string(points1)
      names(spdf) <- "CellID" # Name the column
      
      # Aggregate the points
      spdf@data$points1 <- poly.counts(points1, spdf)
      spdf@data$points2 <- poly.counts(points2, spdf)
      
      # Drop cells with 0 for both counts?
      if (ignore.zeros) {
        spdf <- spdf[which(spdf@data$points1>0 | spdf@data$points2>0),]
        stopifnot( length(which(spdf@data$points1==0 & spdf@data$points2==0)) == 0 )
      }
      
      # Calculate percentages of points in each area (might be useful)
      spdf@data$p1.pct <- 100 * spdf@data$points1 / sum(spdf@data$points1 )
      spdf@data$p2.pct <- 100 * spdf@data$points2 / sum(spdf@data$points2 )
      
      # Calculate the errors 
      
      # Difference in the number of points
      spdf@data$diff <- spdf@data$points1 - spdf@data$points2
      
      # Absolute difference
      spdf@data$abs.diff <- abs(spdf@data$points1 - spdf@data$points2)
      
      # Absolute Difference in percentages
      spdf@data$abs.pct.diff <- abs(spdf@data$p1.pct - spdf@data$p2.pct)
      
      # The Local S Index (slightly more convoluted)
      s <- sppt(base_points.sp = points1, test_points.sp = points2, uoa.sp = spdf) # Calculate the index
      
      # Sanity check - check the sppt package calculates the same percentages as this code
      stopifnot( identical(s$CELLID,spdf$CELLID) ) # Check the cells are in the same order (avoids having to merge on cell ID)
      # XXXX For some reason the two methods aggregat the points very slightly differently. 
      # Uncomment these below, debug, and do something like 'cbind(spdf$points1, s$NumBsePts)' to see the different counts
      #stopifnot(identical(spdf$points1+spdf$points1, s$SumBseTstPts))
      #stopifnot(identical(spdf@data$p1.pct, s$PctBsePts))
      #stopifnot(identical(spdf@data$p2.pct, s$PctTstPts))
      
      # Useful Stats. associated with the S Index 
      # (the whole S object is also returned later too, but these are more convenient to have direct access to)
      spdf@data$localS            <- s@data$localS
      spdf@data$localS.robust     <- s@data$localS.robust
      spdf@data$similarity.robust <- s@data$similarity.robust
      spdf@data$ConfLowP          <- s@data$ConfLowP
      spdf@data$ConfUppP          <- s@data$ConfUppP
      spdf@data$similarity        <- s@data$similarity
      
      # Store this result
      # Formula below (from Nikee) wont work because we might step over some iterations in the inner loop
      #grid.num <- ((i-1) * n.shifts) + j # Calculate total number of grids created so far
      #results[[ grid.num ]] <- spdf
      results[[ counter ]] <- spdf
      counter <- counter + 1
      
      # Now calculte the global errors
     
      # RSS 
      rss <- c(rss, sum( ( spdf@data$points1 - spdf@data$points2 )**2 ))
      # R squared
      r.squared <- c(r.squared, summary(lm(spdf@data$points1 ~ spdf@data$points2, data=spdf@data))$r.squared )
      # RMSE
      rmse <- c(rmse, rmse(spdf@data$points1, spdf@data$points2) )
      # Global S Index (normal and robust). In globalS, each area has same value for global S, so take 1st row
      # arbitrarily. For robust version, need to find the first row that isn't NA (hence use min()).
      #. For
      globalS <- c(globalS, s@data[1,"globalS"]        ) 
      globalS.robust <- c(globalS.robust, s@data[min(which(!is.na(s@data$globalS.robust))),"globalS.robust"] )
      
      # Sometimes useful to keep a reference to the raw results returned by the sppt call (mostly for debugging)
      s.object <- c(s.object, s)
      
      # Also useful to know which iteration number and grid shift this is (useful for naming grids)
      iteration <- c(iteration, i)
      shift <- c(shift, j)
      
      
    } # for shifting grids
    
  } # for cell sizes
  
  
  
  # TODO make parallel: I would need a function that used the outer loop index (i.e. each core will run each resolution)
  # and returned a named list, with each of the components below included in the list. 
  # Then the 'return the results' bit below can go through the results list, exrracting each component
  # (an lapply) and combining them into single vectors
  
  
  
  # Delete the results that used one single large cell as these don't mean anything
  results[[1]] <- NULL
  num.cells <-      num.cells     [2:length(num.cells)]
  cell.areas <-     cell.areas    [2:length(cell.areas)]
  rss <-            rss           [2:length(rss)]
  r.squared <-      r.squared     [2:length(r.squared)]
  rmse <-           rmse          [2:length(rmse)]
  globalS <-        globalS       [2:length(globalS)]
  globalS.robust <- globalS.robust[2:length(globalS.robust)]
  s.object <-       s.object      [2:length(s.object)]
  iteration <-      iteration     [2:length(iteration)]
  shift <-          shift         [2:length(shift)]
  
  # Sanity check - global errors and other info should be vectors of the same length
  stopifnot(
    length(num.cells) == length(cell.areas) &
    length(num.cells) == length(rss) &
    length(num.cells) == length(r.squared) &
    length(num.cells) == length(rmse) & 
    length(num.cells) == length(globalS) & 
    length(num.cells) == length(globalS.robust) &
    length(num.cells) == length(s.object) &
    length(num.cells) == length(iteration) &
    length(num.cells) == length(shift)
  )

  # Return the results
  r <- list(
    "results" = if (return.grids) results else NA,
    "cell.areas" =cell.areas,
    "num.cells" = num.cells,
    "rss" = rss,
    "r.squared" = r.squared,
    "rmse" = rmse,
    "globalS" = globalS,
    "globalS.robust" = globalS.robust,
    "s.object" = if (return.sobject) s.object else NA,
    "iteration" = iteration,
    "shift" = shift
  )
  return(r)
  
} # function


p1 <- BNER[BNER$YEAR==2015,]
p2 <- BNER[BNER$YEAR==2016,]
test.result <- msea(points1 = p1, points2 = p2,
                    N=500, n.shifts=1, ignore.zeros = F, step=10,
                    return.sobject=TRUE,
                    return.grids=TRUE)

save(test.result, file="large_iteration_test.RData")











