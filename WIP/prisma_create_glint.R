#' @title prisma_create_glint
#' @description helper function used to process and save the GLINT data cube
#' @param f input data he5 from caller
#' @param out_file_glint output file name for glint
#' @inheritParams convert_prisma
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
prisma_create_glint <- function(f,
                                out_file_glint,
                                out_format,
                                base_georef){

    message(" - Accessing GLINT dataset - ")

    # Get geo info ----
    geo <- prisma_get_geoloc(f, "1", "HCO")

    glint_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SunGlint_Mask"]][,]
    if (base_georef) {
        rast_glint <- raster::raster(glint_cube, crs = "+proj=longlat +datum=WGS84")
    } else {
        rast_glint <- raster::raster(glint_cube)
    }
    rast_glint <- raster::flip(rast_glint, 1)
    rast_glint[rast_glint == 255] <- NA
    rm(glint_cube)
    gc()

    ex   <- matrix(c(min(geo$lon), max(geo$lon),
                     min(geo$lat), max(geo$lat)),
                   nrow = 2, ncol = 2, byrow = T)
    ex   <- raster::extent(ex)
    if (base_georef) {
        rast_glint <- raster::setExtent(rast_glint, ex, keepres = FALSE)
    }

    message("- Writing GLINT raster -")
    rastwrite_lines(rast_glint, out_file_glint, out_format)
}
