#' @title prisma_create_pan
#' @description helper function used to process and save the PAN data cube
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @param out_file_pan output file name for PAN
#' @inheritParams convert_prisma
#' @return The function is called for its side effects
#' @importFrom hdf5r h5attr
#' @importFrom raster raster extent flip setExtent t
#'
prisma_create_pan <- function(f,
                              proc_lev,
                              source,
                              out_file_pan,
                              out_format,
                              base_georef,
                              fill_gaps,
                              fix_geo){

    # Get geo info ----

    message(" - Accessing PAN dataset - ")
    if (proc_lev %in% c("1")) {
        pan_scale  <- hdf5r::h5attr(f, "ScaleFactor_Pan")
        pan_offset <- hdf5r::h5attr(f, "Offset_Pan")
        pan_cube <- f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source), "/Data Fields/Cube")]][,]
        pan_lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                               "/Geolocation Fields/Latitude")]][,])
        pan_lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                               "/Geolocation Fields/Longitude")]][,])

    } else {
        pan_cube  <- f[[paste0("//HDFEOS/SWATHS/PRS_L", proc_lev, "_PCO/Data Fields/Cube")]][,]
        panscale_min <- hdf5r::h5attr(f, "L2ScalePanMin")
        panscale_max <- hdf5r::h5attr(f, "L2ScalePanMax")
        if (proc_lev == "2D") {
            proj_code <- hdf5r::h5attr(f, "Projection_Id")
            proj_name <- hdf5r::h5attr(f, "Projection_Name")
            xmin  <- hdf5r::h5attr(f, "Product_ULcorner_easting")
            xmax  <- hdf5r::h5attr(f, "Product_LRcorner_easting")
            ymin  <- hdf5r::h5attr(f, "Product_LRcorner_northing")
            ymax  <- hdf5r::h5attr(f, "Product_ULcorner_northing")
            geo <- list(xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax,
                        proj_code = proj_code,
                        proj_name = proj_name)
        }
        if (proc_lev  %in% c("2B", "2C")) {
            pan_lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", gsub("H", "P", source),
                                   "/Geolocation Fields/Latitude")]][,])
            pan_lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", gsub("H", "P", source),
                                   "/Geolocation Fields/Longitude")]][,])
        }
    }

    if (proc_lev %in% c("1", "2B", "2C")) {

        if (base_georef) {
            rast_pan <- raster::raster(pan_cube, crs = "+proj=longlat +datum=WGS84")
            if (proc_lev == "1") {
                rast_pan <- (rast_pan / pan_scale) - pan_offset
            }
            band <- prisma_basegeo(rast_pan, pan_lon, pan_lat, fill_gaps)
        } else {
            rast_pan <- raster::raster(pan_cube)
            rast_pan <- raster::flip(rast_pan, 1)
        }
    } else {
        if (proc_lev == "2D") {
            rast_pan <- raster::raster(pan_cube,
                                       crs = paste0("+proj=utm +zone=", proj_code,
                                                    " +datum=WGS84 +units=m +no_defs"))
            rast_pan <- raster::t(rast_pan)
            ex <- matrix(c(geo$xmin, geo$xmax,
                           geo$ymin, geo$ymax),
                         nrow = 2, ncol = 2, byrow = T)
            ex <- raster::extent(ex)
            if (fix_geo) {
                ex <- ex - 90
            }
            rast_pan <- raster::setExtent(rast_pan, ex, keepres = FALSE)
        }
    }

    rm(pan_cube)
    gc()


    message("- Writing PAN raster -")

    rastwrite_lines(rast_pan, out_file_pan, out_format, proc_lev,
                    scale_min = panscale_min,
                    scale_max = panscale_max)
    rm(rast_pan)
    rm(pan_lon)
    rm(pan_lat)
    gc()
}
