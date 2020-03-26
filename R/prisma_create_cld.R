#' @title prisma_create_cld
#' @description helper function used to process and save the CLOUD data cube
#' @param f input data he5 from caller
#' @param out_file_cld output file name for glint
#' @inheritParams convert_prisma
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
prisma_create_cld <- function(f,
                              out_file_cld,
                              out_format,
                              base_georef){

    message(" - Accessing CLOUD raster - ")

    # Get geo info ----
    geo <- prisma_get_geoloc(f, "1", "HCO")

    cld_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/Cloud_Mask"]][,]
    if (base_georef) {
        rast_cld <- raster::raster(cld_cube, crs = "+proj=longlat +datum=WGS84")
    }
    rast_cld <- raster::flip(rast_cld, 1)
    rast_cld[rast_cld == 255] <- NA
    rm(cld_cube)
    gc()

    ex   <- matrix(c(min(geo$lon), max(geo$lon),
                     min(geo$lat), max(geo$lat)),
                   nrow = 2, ncol = 2, byrow = T)
    ex   <- raster::extent(ex)
    if (base_georef) {
        rast_cld <- raster::setExtent(rast_cld, ex, keepres=F)
    }
    message("- Writing CLD raster -")
    rastwrite_lines(rast_cld, out_file_cld, out_format)
}
