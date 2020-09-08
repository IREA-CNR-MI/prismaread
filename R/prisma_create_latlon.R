#' @title prisma_create_additional
#' @description helper function used to process and save LAT LON datasets
#' @param f input data he5 from caller
#' @param out_file output file name for the dataset
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @inheritParams pr_convert
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
prisma_create_latlon <- function(f,
                                 proc_lev,
                                 out_file,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 in_L2_file = NULL){

    message(" - Accessing LatLon dataset - ")

    # Get geo info ----
    geo <- prisma_get_geoloc(f, proc_lev, "HCO", "VNIR", in_L2_file)

    if (proc_lev != "2D") {
        rast_lat  <- raster::t(raster::raster(geo$lat))
        rast_lon  <- raster::t(raster::raster(geo$lon))
        if (base_georef) {
            rast_lat  <- prisma_basegeo(rast_lat, geo$lon, geo$lat, fill_gaps)
            rast_lon  <- prisma_basegeo(rast_lon, geo$lon, geo$lat, fill_gaps)

        } else {
            rast_lat  <- raster::flip(rast_lat, 1)
            rast_lon  <- raster::flip(rast_lon, 1)
            raster::projection(rast_lat) <- NA
            raster::projection(rast_lon) <- NA
        }
    } else {
        rast_lat  <- raster::raster(
            geo$lat,
            crs = paste0(
                "+proj=utm +zone=", geo$proj_code,
                ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
                " +datum=WGS84 +units=m +no_defs"))
        rast_lon  <- raster::raster(
            geo$lon,
            crs = paste0(
                "+proj=utm +zone=", geo$proj_code,
                ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
                " +datum=WGS84 +units=m +no_defs"))
        # rast_lat  <- raster::t(rast_lat)
        # rast_lon  <- raster::t(rast_lon)
        ex <- matrix(c(geo$xmin - 15, geo$xmin - 15 + dim(rast_lat)[2]*30,
                       geo$ymin - 15, geo$ymin - 15 + dim(rast_lat)[1]*30),
                     nrow = 2, ncol = 2, byrow = T)
        ex <- raster::extent(ex)
        rast_lat  <- raster::setExtent(rast_lat, ex, keepres = FALSE)
        rast_lon  <- raster::setExtent(rast_lon, ex, keepres = FALSE)
    }
    rastang <- raster::stack(rast_lat,
                             rast_lon)
    names(rastang) <- c("lat", "lon")
    gc()
    message(" - Writing LATLON raster - ")
    rastwrite_lines(rastang, out_file, out_format)
    if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file), ".hdr")
        cat("band names = {", paste(names(rastang),collapse=","), "}", "\n",
            file=out_hdr, append=TRUE)
    }
}
