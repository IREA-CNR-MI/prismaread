#' @title pr_create_angles
#' @description helper function used to process and save the ANGLES datasets
#' @param f input data he5 from caller
#' @param out_file output file name for glint
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @inheritParams pr_convert
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
pr_create_angles <- function(f,
                                 proc_lev,
                                 out_file,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 in_L2_file = NULL){


    message(" - Accessing ANGLES dataset - ")

    # Get geo info ----
    geo <- pr_get_geoloc(f, proc_lev, "HCO", "VNIR", in_L2_file)
    if (!is.null(in_L2_file)) {
        f <- try(hdf5r::H5File$new(in_L2_file, mode="r+"))
        proc_lev <- hdf5r::h5attr(f, "Processing_Level")
        if (inherits(f, "try-error")){
            stop("Unable to open the input accessory L2 file as a hdf5 file.
                 Verify your inputs. Aborting!")
        }
    }
    # Process L1/2B/2C data ----
    if (proc_lev != "2D") {

        # Get data from HDF, perform GLT georef if needed and create raster
        # bands
        rast_viewzen    <- raster::raster(
            f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                      "_HCO/Geometric Fields/Observing_Angle")]][,])
        rast_relazang  <- raster::raster(
            f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                      "_HCO/Geometric Fields/Rel_Azimuth_Angle")]][,])
        rast_solzenang <- raster::raster(
            f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                      "_HCO/Geometric Fields/Solar_Zenith_Angle")]][,])
        if (base_georef) {
            rast_viewzen   <- pr_basegeo(rast_viewzen, geo$lon, geo$lat,
                                             fill_gaps)
            rast_relazang  <- pr_basegeo(rast_relazang, geo$lon, geo$lat,
                                             fill_gaps)
            rast_solzenang <- pr_basegeo(rast_solzenang, geo$lon, geo$lat,
            )
        } else {
            rast_viewzen   <- raster::flip(rast_viewzen, 1)
            rast_relazang  <- raster::flip(rast_relazang, 1)
            rast_solzenang <- raster::flip(rast_solzenang, 1)
            raster::projection(rast_viewzen) <- NA
            raster::projection(rast_relazang) <- NA
            raster::projection(rast_solzenang) <- NA
        }
    } else {
        # Process L2D data ----
        #
        # Get data from HDF, create raster bands and set extent ----

        outcrs <- paste0(
            "+proj=utm +zone=", geo$proj_code,
            ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
            " +datum=WGS84 +units=m +no_defs")

        rast_viewzen    <- raster::raster(
            f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                      "_HCO/Geometric Fields/Observing_Angle")]][,],
            crs = outcrs)

        rast_relazang  <- raster::raster(
            f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                      "_HCO/Geometric Fields/Rel_Azimuth_Angle")]][,],
            crs = outcrs)

        rast_solzenang <- raster::raster(
            f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev,
                      "_HCO/Geometric Fields/Solar_Zenith_Angle")]][,],
            crs = outcrs)

        rast_viewzen   <- raster::t(rast_viewzen)
        rast_relazang  <- raster::t(rast_relazang)
        rast_solzenang <- raster::t(rast_solzenang)

        rast_viewzen   <- pr_setext_L2D(geo, rast_viewzen)
        rast_relazang  <- pr_setext_L2D(geo, rast_relazang)
        rast_solzenang <- pr_setext_L2D(geo, rast_solzenang)
    }

    rastang <- raster::stack(rast_viewzen,
                             rast_relazang,
                             rast_solzenang)
    rastang[rastang == 0 ] <- NA
    names(rastang) <- c("view_zenang", "relaz_ang", "solzen_ang")
    gc()
    message(" - Writing ANGLES raster - ")
    pr_rastwrite_lines(rastang, out_file, out_format)

    if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file), ".hdr")
        cat("band names = {", paste(names(rastang),collapse=","), "}", "\n",
            file=out_hdr, append=TRUE)
    }
    rm(rastang ,rast_viewzen, rast_relazang, rast_solzenang)
    rm(geo)
}
