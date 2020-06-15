#' @title prisma_get_geoloc
#' @description helper function used to get geolocation info
#'  from PRISMA data on VNIR and SWIR cubes
#' @return `list` containing required info according to `proc_lev`
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @param wvl `character` "VNIR" or "SWIR" - passed by caller
#' @inheritParams convert_prisma
#' @importFrom hdf5r h5attr
#'
prisma_get_geoloc <- function(f, proc_lev, source, wvl = NULL, in_L2_file = NULL) {
    if (proc_lev == "1") {
        if (is.null(in_L2_file)) {
            if (is.null(wvl) | wvl == "VNIR") {
                lat <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_VNIR")]][,])
                lon <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_VNIR")]][,])
            } else {
                lat <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_SWIR")]][,])
                lon <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_SWIR")]][,])
            }
        } else {
            f2 <- try(hdf5r::H5File$new(in_L2_file, mode="r+"))
            if (inherits(f2, "try-error")){
                stop("Unable to open the input accessory L2 file as a hdf5 file. Verify your inputs. Aborting!")
            }
            proc_lev_f2 <- hdf5r::h5attr(f2, "Processing_Level")
            if (proc_lev_f2 == "1") {
                stop("in_L2_file is not a L2 PRISMA file. Aborting!")
            }
            lat <- raster::t(f2[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev_f2, "_", source, "/Geolocation Fields/Latitude")]][,])
            lon <- raster::t(f2[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev_f2, "_", source, "/Geolocation Fields/Longitude")]][,])
        }
        out <- list(lat = lat, lon = lon)
        return(out)
    } else {
        if (proc_lev == "2D") {
            proj_code <- hdf5r::h5attr(f, "Projection_Id")
            proj_name <- hdf5r::h5attr(f, "Projection_Name")
            proj_epsg <- hdf5r::h5attr(f, "Epsg_Code")
            xmin  <- min(hdf5r::h5attr(f, "Product_ULcorner_easting"),
                         hdf5r::h5attr(f, "Product_LLcorner_easting"))
            xmax  <- max(hdf5r::h5attr(f, "Product_LRcorner_easting"),
                         hdf5r::h5attr(f, "Product_URcorner_easting"))
            ymin  <- min(hdf5r::h5attr(f, "Product_LLcorner_northing"),
                         hdf5r::h5attr(f, "Product_LRcorner_northing"))
            ymax  <- max(hdf5r::h5attr(f, "Product_ULcorner_northing"),
                         hdf5r::h5attr(f, "Product_URcorner_northing"))
            lat <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Latitude")]][,])
            lon <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Longitude")]][,])
            out <- list(xmin = xmin,
                        xmax = xmax,
                        ymin = ymin,
                        ymax = ymax,
                        proj_code = proj_code,
                        proj_name = proj_name,
                        proj_epsg = proj_epsg,
                        lat = lat,
                        lon = lon)
            return(out)
        }

        if (proc_lev  %in% c("2B", "2C")) {
            lat <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Latitude")]][,])
            lon <- raster::t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Longitude")]][,])
            out <- list(lat = lat, lon = lon)
            return(out)
        }
    }
}
