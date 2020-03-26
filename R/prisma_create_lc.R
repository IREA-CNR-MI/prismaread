#' @title prisma_create_glint
#' @description helper function used to process and save the LC data cube
#' @param f input data he5 from caller
#' @param out_file_lc output file name for glint
#' @inheritParams convert_prisma
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent

prisma_create_lc <- function(f,
                             out_file_lc,
                             out_format,
                             base_georef){

    message(" - Accessing LAND COVER raster - ")

    # Get geo info ----
    geo <- prisma_get_geoloc(f, "1", "HCO")

    lc_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SunGlint_Mask"]][,]
    if (base_georef) {
        rast_lc <- raster::raster(lc_cube, crs = "+proj=longlat +datum=WGS84")
    } else {
        rast_lc <- raster::raster(lc_cube)
    }
    rast_lc <- raster::flip(rast_lc, 1)
    rast_lc[rast_lc == 255] <- NA
    rm(lc_cube)
    gc()

    ex   <- matrix(c(min(geo$lon), max(geo$lon),
                     min(geo$lat), max(geo$lat)),
                   nrow = 2, ncol = 2, byrow = T)
    ex   <- raster::extent(ex)
    if (base_georef) {
        rast_lc <- raster::setExtent(rast_lc, ex, keepres = FALSE)
    }

    message("- Writing lc raster -")
    rastwrite_lines(rast_lc, out_file_lc, out_format)
}
