#' @title prisma_create_cld
#' @description helper function used to process and save the VNIR data cube
#' @inheritParams convert_prisma
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom hdf5r h5attr
#' @importFrom raster raster flip extent setExtent stack
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#'
prisma_create_cld <- function(f,
                              out_file_cld,
                              out_format,
                              base_georef){

    message(" - Accessing CLOUD raster - ")

    # Get geo info ----
    geo <- prisma_get_geoloc(f, "1", "HCO")

    cld_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/Cloud_Mask"]][,]
    rast_cld <- raster::raster(cld_cube, crs = "+proj=longlat +datum=WGS84")
    rast_cld <- raster::flip(rast_cld, 1)
    rast_cld[rast_cld == 255] <- NA
    rm(cld_cube)
    gc()

    ex   <- matrix(c(min(geo$lon), max(geo$lon),
                     min(geo$lat), max(geo$lat)),
                   nrow = 2, ncol = 2, byrow = T)
    ex   <- raster::extent(ex)
    rast_cld <- raster::setExtent(rast_cld, ex, keepres=F)

    message("- Writing CLD raster -")
    rastwrite_lines(rast_cld, out_file_cld, out_format)
}
