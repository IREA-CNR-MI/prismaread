#' @title prisma_get_geoloc
#' @description helper function used to get geolocation info
#'  from PRISMA data on VNIR and SWIR cubes
#' @return `list` containing required info according to `proc_lev`
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @inheritParams convert_prisma
#' @importFrom hdf5r h5attr
prisma_get_geoloc <- function(f, proc_lev, source) {

    if (proc_lev == "1") {
        lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_SWIR")]][,])
        lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_SWIR")]][,])
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

            # convert to 32632
            # ll_pt <- sf::st_sf(p = "LL",
            #                    geometry =  sf::st_sfc(sf::st_point(c(xmin, ymin))),
            #                    crs = 4326)
            # ur_pt <- sf::st_sf(p = "UR",
            #                    geometry =  sf::st_sfc(sf::st_point(c(xmax, ymax))),
            #                    crs = 4326)

            # ll_pt <- sf::st_transform(ll_pt, 32632)
            # ur_pt <- sf::st_transform(ur_pt, 32632)
            # ymax  <- hdf5r::h5attr(f, "Product_ULcorner_northing"),
            #            hdf5r::h5attr(f, "Product_URcorner_northing"))
            # out <- list(xmin = sf::st_coordinates(ll_pt)[1],
            #             xmax = sf::st_coordinates(ur_pt)[1],
            #             ymin = sf::st_coordinates(ll_pt)[2],
            #             ymax = sf::st_coordinates(ur_pt)[2],
            #             proj_code = proj_code,
            #             proj_name = proj_name)
            out <- list(xmin = xmin,
                        xmax = xmax,
                        ymin = ymin,
                        ymax = ymax,
                        proj_code = proj_code,
                        proj_name = proj_name)
            return(out)
        }

        if (proc_lev  %in% c("2B", "2C")) {
            xmin  <- min(c(hdf5r::h5attr(f, "Product_LLcorner_long"),
                           hdf5r::h5attr(f, "Product_ULcorner_long"),
                           hdf5r::h5attr(f, "Product_URcorner_long"),
                           hdf5r::h5attr(f, "Product_LRcorner_long")))
            xmax  <- max(c(hdf5r::h5attr(f, "Product_LLcorner_long"),
                           hdf5r::h5attr(f, "Product_ULcorner_long"),
                           hdf5r::h5attr(f, "Product_URcorner_long"),
                           hdf5r::h5attr(f, "Product_LRcorner_long")))
            ymin  <- min(c(hdf5r::h5attr(f, "Product_LLcorner_lat"),
                           hdf5r::h5attr(f, "Product_ULcorner_lat"),
                           hdf5r::h5attr(f, "Product_URcorner_lat"),
                           hdf5r::h5attr(f, "Product_LRcorner_lat")))
            ymax  <- max(c(hdf5r::h5attr(f, "Product_LLcorner_lat"),
                           hdf5r::h5attr(f, "Product_ULcorner_lat"),
                           hdf5r::h5attr(f, "Product_URcorner_lat"),
                           hdf5r::h5attr(f, "Product_LRcorner_lat")))
            lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Latitude")]][,])
            lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L", proc_lev, "_", source, "/Geolocation Fields/Longitude")]][,])
            out <- list(lat = lat, lon = lon,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax)
            return(out)
        }
    }
}
