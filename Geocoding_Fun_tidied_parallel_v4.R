# remove all objects from memory. Good practice if you want to run everything
# from beginning to end, to first clear R's memory just to be absolutely sure
rm(list = ls())

setwd('/Users/nick/mapping/projects/sppt')

# necessary librareis below :)

library(GISTools)
library(rgdal)     
# library(xtable)  
library(parallel) 
library(sppt) 
library(scales) # this is to get transparency in plots
# the next package is a wrapper around the standard parallel functions that 
# shows a progress bar that also estimates how much time is remaining.
# It also works when using multiple CPUs (parallel lapply). It's great!
library(pbapply)


# library(dplyr)

# Load shapefiles
ZIPFILENAME <- "./data-raw/vancouver_public/crime_shp_all_years.zip"
zipfile <- unzip(ZIPFILENAME)
all.crime <- readOGR(dsn="./crime_shp_all_years", layer = "crime_shp_all_years")
# Delete the extracted file
unlink(zipfile)
unlink("crime_shp_all_years", recursive = TRUE) # a left over directory
rm(ZIPFILENAME)

d <- all.crime[all.crime$TYPE=="Break and Enter Residential/Other" & all.crime$YEAR==2016,]
#d <- all.crime[all.crime$TYPE=="Break and Enter Residential/Other",]
rm(all.crime)

#das <- readOGR(dsn = "vancouver_data", layer = "Vancouver_CTs_Census_2016_UTMz10")
das <- readOGR(dsn = "data-raw/vancouver_all", layer = "areas")


# Prepare a vector of the proportions of data to sample (e.g. 0.99, 0.98, ... 0.50)
proportion <- rep(seq(from=0.01, to=0.99, by=0.01), each = 10)


######################################################
#
# the actual function: parallel  --- RUN THIS BLOCK
#
######################################################

set.seed(34637)

# Detect number of cores available
ncores <- detectCores()

# Set up parallel environment using ncores
cl <- makeCluster(ncores, methods = TRUE)

# load necessary packages on all cores
clusterEvalQ(cl, {
  library(sppt)
  # library(dplyr)
}) 

# export objects (in 1 core environment) to the nodes of the cluster
clusterExport(cl, "d")
clusterExport(cl, "das")

# added random seed (appropriately for parallel computing) for replicability
clusterSetRNGStream(cl = cl, iseed = 325636) # I chose 325636

# run all the tests
s.objects.par <- pblapply(proportion, function(p){
  
  samp <- sample(1:nrow(d), size=nrow(d)*p) # A sample of random rows
  base <- d
  test <- d[samp,]
  
  s <- summary.sppt(sppt(base, test, das, bootstrap = TRUE))
  # or should the above line read sppt ?
  
  return(s)
}, cl = cl)

# Close parallel environment
stopCluster(cl)

# extract the first and 2nd element of my list of lists
globalS <-        sapply(s.objects.par, `[[`, 1)



######################################################
#
# plot
#
######################################################

# I've changed the code use a local regression (loess) to get a smooth line,
# two separate plots, and add multiple indicator lines to easy identifying
# geocoding accuracy for a given S value.

# Calculate lines of best fit for global S against proportion (using loess)
fit <-        loess( y ~ x, data = data.frame(x = proportion, y = globalS) ) 

# predicted values
preds <-        predict(fit,        data = data.frame(x = proportion, y = globalS),        se = TRUE)

# what are the x-values for globalS (or globalS.robust) being 0.8, 0.9, 0.95, and 0.99 ?
xvals <- data.frame(svalue = c(0.8, 0.9, 0.95, 0.99), xval = NA)
for(i in 1:nrow(xvals)){
  xvals$xval[i] <-        approx(x = fit$fitted, y = fit$x, xout = xvals$svalue[i])$y
}

# Plot the results: run this as an entire block

# plot globalS
plot(x = proportion, y = globalS, main="Global S", 
     xlab="Geocoding accuracy", ylab="Global S", 
     col = scales::alpha("black", .5), cex = 0.5)

# add best fit line and 95% confidence bands
lines(proportion, preds$fit, col="red")
lines(proportion, preds$fit - qt(0.975, preds$df) * preds$se, lty=2)
lines(proportion, preds$fit + qt(0.975, preds$df) * preds$se, lty=2)

# add lines indicating S value of 0.8, 0.9, 0.95, and 0.99
for (i in 1:nrow(xvals)){
  # check if a line should be drawn to begin with
  if(!is.na(xvals$xval[i])){
    # draw horizonal line
    segments(x0 = proportion[1], y0 = xvals$svalue[i], x1 = xvals$xval[i], y1 = xvals$svalue[i], col = "blue")
    # draw vertical line
    segments(x0 = xvals$xval[i], y0 = 0, x1 = xvals$xval[i], y1 = xvals$svalue[i], col = "blue")
  }
}




#### splines

#install.packages("pspline")
library(pspline)


# penalized splines in pspline
temp2 <- pspline::sm.spline(x = proportion, y = globalS)
fit2 <- predict(temp2, xarg = proportion)

# # what are the x-values for globalS (or globalS.robust) being 0.8, 0.9, 0.95, and 0.99 ?
# xvals <- data.frame(svalue = c(0.8, 0.9, 0.95, 0.99), xval = NA)
# for(i in 1:nrow(xvals)){
#   xvals$xval[i] <-        approx(x = fit2, y = proportion, xout = xvals$svalue[i])$y
# }

# Plot the results: run this as an entire block

# plot globalS
plot(x = proportion, y = globalS, main="Penalized splines", 
     xlab="Geocoding accuracy", ylab="Global S", 
     col = scales::alpha("black", .5), cex = 0.5)

lines(temp2, col="red", lwd = 2) # is kinda wiggly
lines(pspline::sm.spline(x = proportion, y = globalS, df=14), col = "blue", lwd = 2) # much smoother

# # add lines indicating S value of 0.8, 0.9, 0.95, and 0.99
# for (i in 1:nrow(xvals)){
#   # check if a line should be drawn to begin with
#   if(!is.na(xvals$xval[i])){
#     # draw horizonal line
#     segments(x0 = proportion[1], y0 = xvals$svalue[i], x1 = xvals$xval[i], y1 = xvals$svalue[i], col = "blue")
#     # draw vertical line
#     segments(x0 = xvals$xval[i], y0 = 0, x1 = xvals$xval[i], y1 = xvals$svalue[i], col = "blue")
#   }
# }






### ChangePoint detection

# install.packages("ShapeChange")

library(ShapeChange)

par(mfrow=c(2,2))

# penalized in ShapeChange
temp3 <- changept(globalS ~ tp(proportion, sh = 1), pnt = TRUE) 
plot(temp3) # 0.95

# non-penalized in ShapeChange
temp4 <- changept(globalS ~ tp(proportion, sh = 1))
plot(temp4) # 0.92


# giving a 95% bootstrap confidence interval of the changepoint
# and using non-parametric bootstrapping i.o. parametric:
# this takes a while
temp3 <- changept(globalS ~ tp(proportion, sh = 1), pnt = TRUE, ci = TRUE, param = FALSE) 
plot(temp3) # 0.7717-0.99

temp4 <- changept(globalS ~ tp(proportion, sh = 1), ci = TRUE, param = FALSE) 
plot(temp4) # 


