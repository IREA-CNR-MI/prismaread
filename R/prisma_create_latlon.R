#' @title prisma_create_additional
#' @description helper function used to process and save LAT LON datasets
#' @param f input data he5 from caller
#' @param out_file output file name for the dataset
#' @inheritParams convert_prisma
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
prisma_create_latlon <- function(f,
                                 proc_lev,
                                 out_file,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 fix_geo){

    message(" - Accessing LatLon dataset - ")


    # Get geo info ----
    geo <- prisma_get_geoloc(f, proc_lev, "HCO", "VNIR")
    browser()
    rast_lat  <- raster::raster(geo$lat)
    rast_lon  <- raster::raster(geo$lon)

    if (proc_lev != "2D") {
        rast_lat  <- raster::raster(geo$lat)
        rast_lon  <- raster::raster(geo$lon)
        if (base_georef) {
            rast_lat  <- prisma_basegeo(rast_lat, geo$lon, geo$lat, fill_gaps)
            rast_lon  <- prisma_basegeo(rast_lon, geo$lon, geo$lat, fill_gaps)
        } else {
            rast_lat  <- raster::flip(rast_lat, 1)
            rast_lon  <- raster::flip(rast_lon, 1)
        }
    } else {
        rast_lat  <- raster::raster(geo$lat, crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                          " +datum=WGS84 +units=m +no_defs"))
        rast_lon  <- raster::raster(geo$lon, crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                          " +datum=WGS84 +units=m +no_defs"))
        rast_lat  <- raster::t(rast_lat)
        rast_lon  <- raster::t(rast_lon)
        ex <- matrix(c(geo$xmin, geo$xmax,
                       geo$ymin, geo$ymax),
                     nrow = 2, ncol = 2, byrow = T)
        ex <- raster::extent(ex)
        if (fix_geo) {
            ex <- ex - 90
        }
        rast_lat  <- raster::setExtent(rast_lat, ex, keepres = FALSE)
        rast_lon  <- raster::setExtent(rast_lon, ex, keepres = FALSE)
    }
    rastang <- raster::stack(rast_lat,
                             rast_lon)
    names(rastang) <- c("lat", "lon")
    gc()
    message(" - Writing ANGLES raster - ")
    rastwrite_lines(rastang, out_file, out_format)
}
