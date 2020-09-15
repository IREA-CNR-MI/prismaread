#' @title pr_create_additional
#' @description helper function used to process and save LAT LON datasets
#' @param f input data he5 from caller
#' @param out_file output file name for the dataset
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by
#'  caller
#' @inheritParams pr_convert
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
pr_create_latlon <- function(f,
                                 proc_lev,
                                 out_file,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 in_L2_file = NULL){

    message(" - Accessing LatLon dataset - ")

    # Get geo info ----
    geo <- pr_get_geoloc(f, proc_lev, "HCO", "VNIR", in_L2_file)

    if (proc_lev != "2D") {
        rast_lat  <- raster::t(raster::raster(geo$lat))
        rast_lon  <- raster::t(raster::raster(geo$lon))
        if (base_georef) {
            rast_lat  <- pr_basegeo(rast_lat, geo$lon, geo$lat, fill_gaps)
            rast_lon  <- pr_basegeo(rast_lon, geo$lon, geo$lat, fill_gaps)

        } else {
            rast_lat  <- raster::flip(rast_lat, 1)
            rast_lon  <- raster::flip(rast_lon, 1)
            raster::projection(rast_lat) <- NA
            raster::projection(rast_lon) <- NA
        }
    } else {

        outcrs <- paste0(
            "+proj=utm +zone=", geo$proj_code,
            ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
            " +datum=WGS84 +units=m +no_defs")
        rast_lat  <- raster::raster(geo$lat, crs = outcrs)
        rast_lon  <- raster::raster(geo$lon, crs = outcrs)

        rast_lat <- pr_setext_L2D(geo, rast_lat)
        rast_lon <- pr_setext_L2D(geo, rast_lon)

    }

    rastlatlon <- raster::stack(rast_lat, rast_lon)
    names(rastlatlon) <- c("lat", "lon")
    gc()
    message(" - Writing LATLON raster - ")
    pr_rastwrite_lines(rastlatlon, out_file, out_format)

    if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file), ".hdr")
        cat("band names = {", paste(names(rastlatlon),collapse=","), "}", "\n",
            file=out_hdr, append=TRUE)
    }

    rm(rastlatlon ,rast_lat, rast_lon)
    rm(geo)
}
