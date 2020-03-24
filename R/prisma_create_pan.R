#' @title prisma_create_pan
#' @description helper function used to process and save the VNIR data cube
#' @inheritParams convert_prisma
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom hdf5r h5attr
#' @importFrom raster raster flip extent setExtent stack
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#'
prisma_create_pan <- function(f,
                              proc_lev,
                              source,
                              out_file_pan,
                              out_format,
                              base_georef){

    # Get geo info ----

    message(" - Accessing PAN raster - ")
    if (proc_lev == "1") {
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

    if (proc_lev == "1") {

        if (base_georef) {
            rast_pan <- raster::raster(pan_cube, crs = "+proj=longlat +datum=WGS84")
            ex   <- matrix(c(min(pan_lon), max(pan_lon),
                             min(pan_lat), max(pan_lat)),
                           nrow = 2, ncol = 2, byrow = T)
            ex   <- raster::extent(ex)
        } else {
            rast_pan <- raster::raster(pan_cube)
        }
        rast_pan <- raster::flip(rast_pan, 1)

    } else {
        if (proc_lev == "2D") {

            rast_pan <- raster::raster(pan_cube,
                                       crs = paste0("+proj=utm +zone=", proj_code,
                                                    " +datum=WGS84 +units=m +no_defs"))
            rast_pan <- t(rast_pan)
            ex <- matrix(c(geo$xmin, geo$xmax,
                           geo$ymin, geo$ymax),
                         nrow = 2, ncol = 2, byrow = T)
            ex <- raster::extent(ex)
        }
        if (proc_lev %in% c("2B", "2C")) {

            if (base_georef) {
                rast_pan <- raster::raster(pan_cube,
                                           crs = "+proj=longlat +datum=WGS84")
                ex   <- matrix(c(min(pan_lon), max(pan_lon),
                                 min(pan_lat), max(pan_lat)),
                               nrow = 2, ncol = 2, byrow = T)
                ex <- raster::extent(ex)

            } else {
                rast_pan <- raster::raster(pan_cube)
            }
            rast_pan <- raster::flip(rast_pan, 1)
        }
    }

    rm(pan_cube)
    gc()

    if (base_georef | proc_lev == "2D") {
        rast_pan <- raster::setExtent(rast_pan, ex, keepres = FALSE)
    }

    message("- Writing PAN raster -")

    rastwrite_lines(rast_pan, out_file_pan, out_format, proc_lev,
                    scale_min = panscale_min,
                    scale_max = panscale_max)
    rm(rast_pan)
    rm(pan_lon)
    rm(pan_lat)
    gc()
}
