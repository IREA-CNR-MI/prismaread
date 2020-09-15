#' @title pr_create_pan
#' @description helper function used to process and save the PAN data cube
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by
#'  caller
#' @param out_file_pan output file name for PAN
#' @inheritParams pr_convert
#' @return The function is called for its side effects
#' @importFrom hdf5r h5attr
#' @importFrom raster raster extent flip setExtent t
#'
pr_create_pan <- function(f,
                          proc_lev,
                          source,
                          out_file_pan,
                          out_format,
                          base_georef,
                          fill_gaps,
                          in_L2_file = NULL){

    # Get geo info ----
    geo <- pr_get_geoloc(f, proc_lev, source, wvl = "PAN", in_L2_file)

    message(" - Accessing PAN dataset - ")
    if (proc_lev %in% c("1")) {
        pan_scale  <- hdf5r::h5attr(f, "ScaleFactor_Pan")
        pan_offset <- hdf5r::h5attr(f, "Offset_Pan")
        pan_cube   <- f[[paste0("/HDFEOS/SWATHS/PRS_L1_",
                                gsub("H", "P", source),
                                "/Data Fields/Cube")]][,]
    } else {
        pan_cube  <- f[[paste0("//HDFEOS/SWATHS/PRS_L", proc_lev,
                               "_PCO/Data Fields/Cube")]][,]
        panscale_min <- hdf5r::h5attr(f, "L2ScalePanMin")
        panscale_max <- hdf5r::h5attr(f, "L2ScalePanMax")
    }

    if (proc_lev %in% c("1", "2B", "2C")) {

        if (base_georef) {
            message("Applying bowtie georeferencing")
            rast_pan <- raster::raster(pan_cube,
                                       crs = "+proj=longlat +datum=WGS84")
            if (proc_lev == "1") {
                rast_pan <- (rast_pan / pan_scale) - pan_offset
            }
            rast_pan <- pr_basegeo(rast_pan, geo$lon, geo$lat, fill_gaps)
        } else {
            rast_pan <- raster::raster(pan_cube)
            rast_pan <- raster::flip(rast_pan, 1)
            raster::projection(rast_pan) <- NA
        }
    } else {
        if (proc_lev == "2D") {
            rast_pan <- raster::raster(
                pan_cube,
                crs = paste0("+proj=utm +zone=", geo$proj_code,
                             ifelse(substring(geo$proj_epsg, 3, 3) == 7,
                                    " +south", ""),
                             " +datum=WGS84 +units=m +no_defs"))
            rast_pan <- raster::t(rast_pan)
            ex <- matrix(c(geo$xmin - 2.5, geo$xmin - 2.5 + dim(rast_pan)[2]*5,
                           geo$ymin - 2.5, geo$ymin - 2.5 + dim(rast_pan)[1]*5),
                         nrow = 2, ncol = 2, byrow = TRUE)
            ex <- raster::extent(ex)
            rast_pan <- raster::setExtent(rast_pan, ex, keepres = FALSE)
        }
    }

    rm(pan_cube)
    gc()

    message("- Writing PAN raster -")

    pr_rastwrite_lines(rast_pan, out_file_pan, out_format, proc_lev,
                       scale_min = panscale_min,
                       scale_max = panscale_max)
    rm(rast_pan)
    rm(geo)
    gc()
}
