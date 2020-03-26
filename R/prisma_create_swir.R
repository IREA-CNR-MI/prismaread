#' @title prisma_create_swir
#' @description helper function used to process and save the SWIR data cube
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @param out_file_swir output file name for SWIR
#' @param wl_swir passed by caller - array of PRISMA SWIR wavelengths
#' @param order_swir passed by caller - ordering of array of PRISMA SWIR wavelengths
#' @param fwhm_swir passed by caller - array of PRISMA SWIR fwhms
#' @inheritParams convert_prisma
#' @return the function is called for its side effects
#' @importFrom hdf5r h5attr
#' @importFrom raster raster extent flip t setExtent stack
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#'
prisma_create_swir <- function(f,
                               proc_lev,
                               source,
                               out_file_swir,
                               out_format,
                               base_georef,
                               wl_swir,
                               order_swir,
                               fwhm_swir,
                               apply_errmatrix,
                               ERR_MATRIX){

    # Get geo info ----
    geo <- prisma_get_geoloc(f, proc_lev, source)

    if (proc_lev == "1") {
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source,
                               "/Data Fields/SWIR_Cube")]][,,]
        swir_scale  <- hdf5r::h5attr(f, "ScaleFactor_Swir")
        swir_offset <- hdf5r::h5attr(f, "Offset_Swir")
    } else {
        swir_max <- hdf5r::h5attr(f, "L2ScaleSwirMax")
        swir_min <- hdf5r::h5attr(f, "L2ScaleSwirMin")
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                               source, "/Data Fields/SWIR_Cube")]][,,]
    }

    ind_band <- 1
    for (band_swir in 1:173) {
        if (wl_swir[band_swir] != 0) {
            if(proc_lev == "1") {
                if (base_georef) {
                    band <- raster::raster((swir_cube[,order_swir[band_swir], ]),
                                           crs = "+proj=longlat +datum=WGS84")
                    ex <- matrix(c(min(geo$lon), max(geo$lon),
                                   min(geo$lat), max(geo$lat)),
                                 nrow = 2, ncol = 2, byrow = T)
                    ex <- raster::extent(ex)
                } else {
                    band <- raster::raster((swir_cube[,order_swir[band_swir], ]))
                }
                band <- (band / swir_scale) - swir_offset
                band <- raster::flip(band, 1)
            } else {

                if (proc_lev == "2D") {
                    band <- raster::raster((swir_cube[,order_swir[band_swir], ]),
                                           crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                        " +datum=WGS84 +units=m +no_defs"))
                    band <- raster::t(band)
                    ex <- matrix(c(geo$xmin, geo$xmax,
                                   geo$ymin, geo$ymax),
                                 nrow = 2, ncol = 2, byrow = T)
                    ex <- raster::extent(ex)
                }

                if (proc_lev %in% c("2B", "2C")) {
                    if (base_georef) {
                        band <- raster::raster((swir_cube[,order_swir[band_swir], ]),
                                               crs = "+proj=longlat +datum=WGS84")
                        ex <- matrix(c(min(geo$lon), max(geo$lon),
                                       min(geo$lat), max(geo$lat)),
                                     nrow = 2, ncol = 2, byrow = T)
                        ex <- raster::extent(ex)
                    } else {
                        band <- raster::raster(swir_cube[,order_swir[band_swir], ])
                    }
                    band <- raster::flip(band, 1)

                }
            }
            if (base_georef | proc_lev == "2D") {
                band <- raster::setExtent(band, ex, keepres = FALSE)
            }

            if (ind_band == 1) {
                rast_swir <- band
            } else {
                rast_swir <- raster::stack(rast_swir, band)
            }
            ind_band <- ind_band + 1
        }
    }

    wl_swir   <- wl_swir[wl_swir != 0]
    fwhm_swir <- fwhm_swir[fwhm_swir != 0]
    names(rast_swir) <- paste0("wl_", round(wl_swir, digits = 4))
    rm(swir_cube)
    rm(band)
    gc()
    message("- Writing SWIR raster -")
    rastwrite_lines(rast_swir,
                    out_file_swir,
                    out_format,
                    proc_lev,
                    scale_min = swir_min,
                    scale_max = swir_max)
    if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file_swir), ".hdr")

        write(c("wavelength = {",
                paste(round(wl_swir, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write(c("fwhm = {",
                paste(round(fwhm_swir, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
    }

    out_file_txt <- paste0(tools::file_path_sans_ext(out_file_swir), "_meta.txt")
    utils::write.table(data.frame(band = 1:length(wl_swir),
                                  wl   = wl_swir,
                                  fwhm = fwhm_swir, stringsAsFactors = FALSE),
                       file = out_file_txt, row.names = FALSE)

}
