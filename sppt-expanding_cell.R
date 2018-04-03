#' @import sppt
NULL

#' An extension to the sppt package that runs the test over a number of grids with varying
#' resolutions.
#' TODO - integrate this script into the sppt package properly.
#' @param min.area The minimum square area (in square meters) of the finest (highest resolution)
#' grid.
#' @param num.scales The number of times that the resolution should be increased. The smallest
#' resolution will be equal to `min.area`, and the largest is one cell over the entire area.
sppt.exp.cell <- function(base_points.sp, test_points.sp, uoa.sp, 
                          min.area = 100, num.scales = 20,
                          outputlist=FALSE, nsamples=200, percpoints=85, conf_level=95){
    
}
                          
