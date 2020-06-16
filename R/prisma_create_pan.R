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
                              fix_geo,
                              in_L2_file = NULL){

    # Get geo info ----

    message(" - Accessing PAN dataset - ")
    if (proc_lev %in% c("1")) {
        pan_scale  <- hdf5r::h5attr(f, "ScaleFactor_Pan")
        pan_offset <- hdf5r::h5attr(f, "Offset_Pan")
        pan_cube <- f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source), "/Data Fields/Cube")]][,]
        if (is.null(in_L2_file)){
            pan_lat <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                                           "/Geolocation Fields/Latitude")]][,])
            pan_lon <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                                           "/Geolocation Fields/Longitude")]][,])
        } else {
            f2 <- try(hdf5r::H5File$new(in_L2_file, mode="r+"))
            if (inherits(f2, "try-error")){
                stop("Unable to open the input accessory L2 file as a hdf5 file. Verify your inputs. Aborting!")
            }
            proc_lev_f2 <- hdf5r::h5attr(f2, "Processing_Level")
            if (proc_lev_f2 == "1") {
                stop("in_L2_file is not a L2 PRISMA file. Aborting!")
            }
            pan_lat <- raster::t(f2[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev_f2, "_", gsub("H", "P", source),
                                            "/Geolocation Fields/Latitude")]][,])
            pan_lon <- raster::t(f2[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev_f2, "_", gsub("H", "P", source),
                                            "/Geolocation Fields/Longitude")]][,])
        }

    } else {
        pan_cube  <- f[[paste0("//HDFEOS/SWATHS/PRS_L", proc_lev, "_PCO/Data Fields/Cube")]][,]
        panscale_min <- hdf5r::h5attr(f, "L2ScalePanMin")
        panscale_max <- hdf5r::h5attr(f, "L2ScalePanMax")
        if (proc_lev == "2D") {
            proj_code <- hdf5r::h5attr(f, "Projection_Id")
            proj_name <- hdf5r::h5attr(f, "Projection_Name")
            proj_epsg <- hdf5r::h5attr(f, "Epsg_Code")
            xmin  <- hdf5r::h5attr(f, "Product_ULcorner_easting")
            xmax  <- hdf5r::h5attr(f, "Product_LRcorner_easting")
            ymin  <- hdf5r::h5attr(f, "Product_LRcorner_northing")
            ymax  <- hdf5r::h5attr(f, "Product_ULcorner_northing")
            geo <- list(xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax,
                        proj_code = proj_code,
                        proj_name = proj_name,
                        proj_epsg = proj_epsg)
        }
        if (proc_lev  %in% c("2B", "2C")) {
            pan_lat <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", gsub("H", "P", source),
                                           "/Geolocation Fields/Latitude")]][,])
            pan_lon <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", gsub("H", "P", source),
                                           "/Geolocation Fields/Longitude")]][,])
        }
    }

    if (proc_lev %in% c("1", "2B", "2C")) {

        if (base_georef) {
            message("Applying bowtie georeferencing")
            rast_pan <- raster::raster(pan_cube, crs = "+proj=longlat +datum=WGS84")
            if (proc_lev == "1") {
                rast_pan <- (rast_pan / pan_scale) - pan_offset
            }
            rast_pan <- prisma_basegeo(rast_pan, pan_lon, pan_lat, fill_gaps)
        } else {
            rast_pan <- raster::raster(pan_cube)
            rast_pan <- raster::flip(rast_pan, 1)
        }
    } else {
        if (proc_lev == "2D") {
            rast_pan <- raster::raster(pan_cube,
                                       crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                    ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
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
