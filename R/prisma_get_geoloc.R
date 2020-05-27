#' @title prisma_get_geoloc
#' @description helper function used to get geolocation info
#'  from PRISMA data on VNIR and SWIR cubes
#' @return `list` containing required info according to `proc_lev`
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @inheritParams convert_prisma
#' @importFrom hdf5r h5attr
prisma_get_geoloc <- function(f, proc_lev, source, wvl = NULL) {

    if (proc_lev == "1") {
        if (is.null(wvl) | wvl == "VNIR") {
            lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_VNIR")]][,])
            lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_VNIR")]][,])
        } else {
            lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_SWIR")]][,])
            lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_SWIR")]][,])
        }
        out <- list(lat = lat, lon = lon)
        return(out)
    } else {
        if (proc_lev == "2D") {

            proj_code <- hdf5r::h5attr(f, "Projection_Id")
            proj_name <- hdf5r::h5attr(f, "Projection_Name")
            xmin  <- min(hdf5r::h5attr(f, "Product_ULcorner_easting"),
                         hdf5r::h5attr(f, "Product_LLcorner_easting"))
            xmax  <- max(hdf5r::h5attr(f, "Product_LRcorner_easting"),
                         hdf5r::h5attr(f, "Product_URcorner_easting"))
            ymin  <- min(hdf5r::h5attr(f, "Product_LLcorner_northing"),
                         hdf5r::h5attr(f, "Product_LRcorner_northing"))
            ymax  <- max(hdf5r::h5attr(f, "Product_ULcorner_northing"),
                         hdf5r::h5attr(f, "Product_URcorner_northing"))
            lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Latitude")]][,])
            lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Longitude")]][,])

            out <- list(xmin = xmin,
                        xmax = xmax,
                        ymin = ymin,
                        ymax = ymax,
                        proj_code = proj_code,
                        proj_name = proj_name,
                        lat = lat,
                        lon = lon)
            return(out)
        }

        if (proc_lev  %in% c("2B", "2C")) {
            lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Latitude")]][,])
            lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Longitude")]][,])
            out <- list(lat = lat, lon = lon)
            return(out)
        }
    }
}
