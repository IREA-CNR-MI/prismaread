#' @title prisma_basegeo
#' @description Function used to create a georeferenced raster from L1/L2B/L2C datasets,
#'  based on the Latitude and Longitude data, and exploiting the procedure described in
#'  https://www.harrisgeospatial.com/docs/backgroundgltbowtiecorrection.html
#' @param band `Raster` layer to be georefewrenced
#' @param lon `Matrix` containing the longitudes of each pixel (as derived from the
#'  geolocation fields datasets of the hdf5)
#' @param lat `Matrix` containing the latitudes of each pixel (as derived from the
#'  geolocation fields datasets of the hdf5)
#' @param fill_gaps `logical` If TRUE, pixels with no values in the georeferenced image are
#'  filled based on a 3x3 average focal filter of the neighbouring valid pixels, Default: TRUE
#' @return `Raster` georeferenced dataset
#' @details The function is based on the "GLT and Bowtie Correction" technique used in ENVI, and
#' described in https://www.harrisgeospatial.com/docs/backgroundgltbowtiecorrection.html.
#' Note that the 7x7 interpolation step for pixels still missing aftere the 3x3 interpolation is not (yet)
#' implemnted, as well as the nearest neighbour interpolation.
#' @seealso
#'  \code{\link[raster]{getValues}}
#' @rdname prisma_basegeo
#' @importFrom raster values raster focal
#' @importFrom stats median
prisma_basegeo <- function(band, lon, lat, fill_gaps = TRUE) {

    # https://www.harrisgeospatial.com/docs/backgroundgltbowtiecorrection.html
    #
    # Estimate the median X and Y pixel sizes, using the center column and row
    # from the GLT.
    numlines <- dim(lon)[1]
    numcols  <- dim(lat)[1]

    psize_x  <- abs(stats::median(diff(lon[, round(numcols / 2) ], lag = 1)))
    psize_y  <- abs(stats::median(diff(lat[round(numlines / 2), ], lag = 1)))

    # Compute the size of the grid, using the following IDL notation.
    # The CEIL function returns the closest integer greater than or equal to
    # its argument.

    ncols <- ceiling((max(lon, na.rm = TRUE) - min(lon, na.rm = TRUE)) / psize_x)
    nrows <- ceiling((max(lat, na.rm = TRUE) - min(lat, na.rm = TRUE)) / psize_y)

    # Set the grid ground control point (GCP) as [0, 0, MinX, MaxY].
    # Map all X and Y entries in the GLT to the output grid, excluding those that
    # are flagged as "bad." If the image crosses the International Date Line, a
    # value of 360 is added to X coordinates that are less than 0. This mapping
    # can be represented as an index into the output grid.
    out_grd <- matrix(data = NA, nrow = nrows, ncol = ncols)
    minlon  <- min(lon, na.rm = TRUE)
    minlat  <- min(lat, na.rm = TRUE)
    vals    <- raster::values(band)

    columns <- round((lon - minlon) / psize_x) + 1
    rows    <- nrows - round((lat - minlat) / psize_y)

    # columns <- round((lon - minlon) / psize_x) + 1
    # rows    <- nrows - round((lat - minlat) / psize_y) - 1

    #♫ remove data if out of bounds to avoid potential crashes
    columns[columns > ncols] <- 0
    rows[rows > nrows] <- 0

    for (indpix in 1:(numlines*numcols)) {
        out_grd[rows[indpix], columns[indpix]] <- vals[indpix]
    }

    outrast <- raster::raster(out_grd, xmn = min(lon, na.rm = TRUE), xmx = min(lon, na.rm = TRUE) + ncols*psize_x,
                      ymn = max(lat, na.rm = TRUE) - nrows*psize_y, ymx = max(lat, na.rm = TRUE),
                      crs = "+init=epsg:4326")

    if (fill_gaps) {
        message("Filling Gaps")
        outrast <- raster::focal(outrast, w=matrix(1,3,3), fun=mean, na.rm = TRUE, NAonly = TRUE)
    }

    return(outrast)

}

