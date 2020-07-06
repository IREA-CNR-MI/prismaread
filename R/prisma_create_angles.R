#' @title prisma_create_angles
#' @description helper function used to process and save the ANGLES datasets
#' @param f input data he5 from caller
#' @param out_file output file name for glint
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @inheritParams convert_prisma
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
prisma_create_angles <- function(f,
                                 proc_lev,
                                 out_file,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 fix_geo,
                                 in_L2_file = NULL){

    message(" - Accessing ANGLES dataset - ")

    # Get geo info ----
    geo <- prisma_get_geoloc(f, proc_lev, "HCO", "VNIR", in_L2_file)
    if (!is.null(in_L2_file)) {
        f <- try(hdf5r::H5File$new(in_L2_file, mode="r+"))
        proc_lev <- hdf5r::h5attr(f, "Processing_Level")
        if (inherits(f, "try-error")){
            stop("Unable to open the input accessory L2 file as a hdf5 file. Verify your inputs. Aborting!")
        }
    }

    if (proc_lev != "2D") {
        rast_obsang    <- raster::raster(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                                                   "_HCO/Geometric Fields/Observing_Angle")]][,])
        rast_relazang  <- raster::raster(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                                                   "_HCO/Geometric Fields/Rel_Azimuth_Angle")]][,])
        rast_solzenang <- raster::raster(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                                                   "_HCO/Geometric Fields/Solar_Zenith_Angle")]][,])
        if (base_georef) {
            rast_obsang    <- prisma_basegeo(rast_obsang, geo$lon, geo$lat, fill_gaps)
            rast_relazang  <- prisma_basegeo(rast_relazang, geo$lon, geo$lat, fill_gaps)
            rast_solzenang <- prisma_basegeo(rast_solzenang, geo$lon, geo$lat, fill_gaps)
        } else {
            rast_obsang    <- raster::flip(rast_obsang, 1)
            rast_relazang  <- raster::flip(rast_relazang, 1)
            rast_solzenang <- raster::flip(rast_solzenang, 1)
        }
    } else {
        rast_obsang    <- raster::raster(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                                                   "_HCO/Geometric Fields/Observing_Angle")]][,],
                                         crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                      ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
                                                      " +datum=WGS84 +units=m +no_defs"))
        rast_relazang  <- raster::raster(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                                                   "_HCO/Geometric Fields/Rel_Azimuth_Angle")]][,],
                                         crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                      ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
                                                      " +datum=WGS84 +units=m +no_defs"))
        rast_solzenang <- raster::raster(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                                                   "_HCO/Geometric Fields/Solar_Zenith_Angle")]][,],
                                         crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                      ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
                                                      " +datum=WGS84 +units=m +no_defs"))
        rast_obsang    <- raster::t(rast_obsang)
        rast_relazang  <- raster::t(rast_relazang)
        rast_solzenang <- raster::t(rast_solzenang)
        ex <- matrix(c(geo$xmin - 15, geo$xmin - 15 + dim(rast_obsang)[2]*30,
                       geo$ymin - 15, geo$ymin - 15 + dim(rast_obsang)[1]*30),
                     nrow = 2, ncol = 2, byrow = T)
        ex <- raster::extent(ex)
        rast_obsang    <- raster::setExtent(rast_obsang, ex, keepres = FALSE)
        rast_relazang  <- raster::setExtent(rast_relazang, ex, keepres = FALSE)
        rast_solzenang <- raster::setExtent(rast_solzenang, ex, keepres = FALSE)

    }
    rastang <- raster::stack(rast_obsang,
                             rast_relazang,
                             rast_solzenang)
    rastang[rastang == 0 ] <- NA
    names(rastang) <- c("obs_ang", "relaz_ang", "solzen_ang")
    gc()
    message(" - Writing ANGLES raster - ")
    rastwrite_lines(rastang, out_file, out_format)
    if (out_format == "ENVI") {
        cat("band names = {", paste(names(rastang),collapse=","), "}", "\n",
            file=raster::extension(out_file, "hdr"), append=TRUE)
    }
}
