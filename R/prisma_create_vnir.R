#' @title prisma_create_vnir
#' @description helper function used to process and save the VNIR data cube
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by caller
#' @param out_file_vnir output file name for VNIR
#' @param wl_vnir passed by caller - array of PRISMA VNIR wavelengths
#' @param order_vnir passed by caller - ordering of array of PRISMA VNIR wavelengths
#' @param fwhm_vnir passed by caller - array of PRISMA VNIR fwhms
#' @inheritParams convert_prisma
#' @return the function is called for its side effects
#' @importFrom hdf5r h5attr
#' @importFrom raster raster extent flip t setExtent stack
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#'
prisma_create_vnir <- function(f,
                               proc_lev,
                               source,
                               out_file_vnir,
                               out_format,
                               base_georef,
                               wl_vnir,
                               order_vnir,
                               fwhm_vnir,
                               apply_errmatrix,
                               ERR_MATRIX){

    # Get geo info ----
    geo <- prisma_get_geoloc(f, proc_lev, source)

    if(proc_lev == 1) {
        vnir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source, "/Data Fields/VNIR_Cube")]][,,]
        vnir_scale  <- hdf5r::h5attr(f, "ScaleFactor_Vnir")
        vnir_offset <- hdf5r::h5attr(f, "Offset_Vnir")
        if (any(apply_errmatrix | ERR_MATRIX)) {
            err_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source,
                                  "/Data Fields/VNIR_PIXEL_SAT_ERR_MATRIX/")]]
        }
    } else {
        vnir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                               source, "/Data Fields/VNIR_Cube")]][,,]
        vnir_max <- hdf5r::h5attr(f, "L2ScaleVnirMax")
        vnir_min <- hdf5r::h5attr(f, "L2ScaleVnirMin")
        if (any(apply_errmatrix | ERR_MATRIX)) {
            err_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                                  source, "/Data Fields/VNIR_PIXEL_SAT_ERR_MATRIX")]][,,]
        }
    }
    ind_vnir <- 1
    for (band_vnir in 1:66) {
        if (wl_vnir[band_vnir] != 0) {

            if(proc_lev == "1") {
                if (base_georef) {
                    band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]),
                                           crs = "+proj=longlat +datum=WGS84")
                    ex <- matrix(c(min(geo$lon), max(geo$lon),
                                   min(geo$lat), max(geo$lat)), nrow = 2, ncol = 2, byrow = T)
                    ex <- raster::extent(ex)
                } else {
                    band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]))
                }
                band <- (band / vnir_scale) - vnir_offset
                band <- raster::flip(band, 1)

            } else {
                if (proc_lev == "2D") {
                    band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]),
                                           crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                        " +datum=WGS84 +units=m +no_defs"))
                    band <- raster::t(band)
                    ex <- matrix(c(geo$xmin, geo$xmax,  geo$ymin, geo$ymax),
                                 nrow = 2, ncol = 2, byrow = T)
                    ex <- raster::extent(ex)
                }
                if (proc_lev %in% c("2B", "2C")) {
                    if (base_georef) {
                        band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]),
                                               crs = "+proj=longlat +datum=WGS84")
                        ex <- matrix(c(min(geo$lon), max(geo$lon),
                                       min(geo$lat), max(geo$lat)),
                                     nrow = 2, ncol = 2, byrow = T)
                        ex <- raster::extent(ex)
                        ex <- ex + 30
                    } else {
                        band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]))
                    }
                    band <- flip(band, 1)
                    }
                }

                if (base_georef | proc_lev == "2D") {
                    band <- raster::setExtent(band, ex, keepres <- FALSE)
                }
                # res(band) <- c(30,30)
                if (ind_vnir == 1) {
                    rast_vnir <- band
                } else {
                    rast_vnir <- raster::stack(rast_vnir, band)
                }
                ind_vnir <- ind_vnir + 1
            }
        }
        # res(rast_vnir) <- c(30,30)
        wl_vnir   <- wl_vnir[wl_vnir != 0]
        fwhm_vnir <- fwhm_vnir[fwhm_vnir != 0]
        names(rast_vnir) <- paste0("wl_", round(wl_vnir, digits = 4))
        rm(vnir_cube)
        rm(band)
        gc()


        message("- Writing VNIR raster -")

        rastwrite_lines(rast_vnir,
                        out_file_vnir,
                        out_format,
                        proc_lev,
                        scale_min = vnir_min,
                        scale_max = vnir_max)
        if (out_format == "ENVI") {

            out_hdr <- paste0(tools::file_path_sans_ext(out_file_vnir), ".hdr")
            write(c("wavelength = {",
                    paste(round(wl_vnir, digits = 4), collapse = ","), "}"),
                  out_hdr, append = TRUE)
            write(c("fwhm = {",
                    paste(round(fwhm_vnir, digits = 4), collapse = ","), "}"),
                  out_hdr, append = TRUE)
        }

        out_file_txt <- paste0(tools::file_path_sans_ext(out_file_vnir), "_meta.txt")
        utils::write.table(data.frame(band = 1:length(wl_vnir),
                                      wl   = wl_vnir,
                                      fwhm = fwhm_vnir, stringsAsFactors = FALSE),
                           file = out_file_txt, row.names = FALSE)
    }
