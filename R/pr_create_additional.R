#' @title pr_create_additional
#' @description helper function used to process and save additional data sets
#'   such as CLOUD, LC and GLINT
#' @param f input data he5 from caller
#' @param type `character` type of dataset to be created ("CLD", "LC" or
#'   "GLINT")
#' @param out_file output file name for the dataset
#' @inheritParams pr_convert
#' @return The function is called for its side effects
#' @importFrom raster raster flip extent setExtent
#'
pr_create_additional <- function(f,
                                     type,
                                     out_file,
                                     out_format,
                                     base_georef,
                                     fill_gaps,
                                     in_L2_file = NULL){

    message(" - Accessing ", type, " dataset - ")

    # Get geo info ----
    geo <- pr_get_geoloc(f, "1", "HCO", "VNIR", in_L2_file)

    cube <- switch(
        type,
        "CLD"  = f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/Cloud_Mask"]][,],
        "LC"   = f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/LandCover_Mask"]][,],
        "GLINT" = f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SunGlint_Mask"]][,],
    )

    rast <- raster::raster(cube)
    if (base_georef) {
        message("Applying bowtie georeferencing")
        rast <- pr_basegeo(rast, geo$lon, geo$lat, fill_gaps)
    } else {
        rast <- raster::flip(rast, 1)
        raster::projection(rast) <- NA
    }

    if (type %in% c("CLD", "LC", "GLINT")) {
        rast[rast == 255] <- NA
    }
    rm(cube)
    gc()
    message(" - Writing  ", type, " raster - ")
    pr_rastwrite_lines(rast, out_file, out_format)
}
