#' @title pr_create_swir
#' @description helper function used to process and save the SWIR data cube
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by
#'  caller
#' @param out_file_swir output file name for SWIR
#' @param wl_swir passed by caller - array of PRISMA SWIR wavelengths
#' @param order_swir passed by caller - ordering of array of PRISMA SWIR
#'  wavelengths
#' @param fwhm_swir passed by caller - array of PRISMA SWIR fwhms
#' @inheritParams pr_convert
#' @return the function is called for its side effects
#' @importFrom hdf5r h5attr
#' @importFrom raster raster extent flip t setExtent stack
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#'
pr_create_swir <- function(f,
                               proc_lev,
                               source,
                               out_file_swir,
                               out_format,
                               base_georef,
                               fill_gaps,
                               wl_swir,
                               order_swir,
                               fwhm_swir,
                               apply_errmatrix,
                               ERR_MATRIX,
                               selbands_swir = NULL,
                               in_L2_file = NULL){

    # Get geo info ----
    geo <- pr_get_geoloc(f, proc_lev, source, wvl = "SWIR", in_L2_file)

    # Get the datacube and required attributes frim hdr ----
    if (proc_lev == "1") {
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source,
                               "/Data Fields/SWIR_Cube")]][,,]
        swir_scale  <- hdf5r::h5attr(f, "ScaleFactor_Swir")
        swir_offset <- hdf5r::h5attr(f, "Offset_Swir")
        if (any(apply_errmatrix | ERR_MATRIX)) {
            err_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source,
                                  "/Data Fields/SWIR_PIXEL_SAT_ERR_MATRIX/")]]
        }
    } else {
        swir_max  <- hdf5r::h5attr(f, "L2ScaleSwirMax")
        swir_min  <- hdf5r::h5attr(f, "L2ScaleSwirMin")
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                               source, "/Data Fields/SWIR_Cube")]][,,]
        if (any(apply_errmatrix | ERR_MATRIX)) {
            err_cube <- f[[paste0(
                "HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                source, "/Data Fields/SWIR_PIXEL_L2_ERR_MATRIX")]][,,]
        }
    }

    # Get the different bands in order of wvl, and convert to `raster`  ----
    # bands. Also georeference if needed
    ind_band <- 1
    if (is.null(selbands_swir)) {
        seqbands <- 1:173
    } else {
        seqbands <- unlist(lapply(
            selbands_swir, FUN = function(x) which.min(abs(wl_swir - x))))
    }

    for (band_swir in seqbands) {

        if (wl_swir[band_swir] != 0) {
            if (proc_lev %in% c("1", "2B", "2C")) {
                if (base_georef) {
                    message("Importing Band: ", band_swir,
                            " (", wl_swir[band_swir], ") of: 173 and
                            applying bowtie georeferencing")
                    band <- raster::raster(
                        (swir_cube[,order_swir[band_swir], ]),
                        crs = "+proj=longlat +datum=WGS84")
                    lat  <- geo$lat
                    lon  <- geo$lon
                    if (proc_lev == "1") {
                        band <- (band / swir_scale) - swir_offset
                    }
                    band <- pr_basegeo(band, lon, lat, fill_gaps)
                    if (ERR_MATRIX | apply_errmatrix) {
                        satband <- raster::raster(
                            (err_cube[,order_swir[band_swir], ]),
                            crs = "+proj=longlat +datum=WGS84")
                        satband <- pr_basegeo(satband, lon, lat, fill_gaps)
                    }
                    if (apply_errmatrix) {
                        band[satband > 0] <- NA
                    }
                } else {
                    message("Importing Band: ", band_swir, " (",
                            wl_swir[band_swir], ") of: of: 173")
                    band <- raster::raster(swir_cube[,order_swir[band_swir], ],
                                           crs = "+proj=longlat +datum=WGS84")
                    if (proc_lev == "1") {
                        band <- (band / swir_scale) - swir_offset
                    }
                    band <- raster::flip(band, 1)
                    raster::projection(band) <- NA

                    if (ERR_MATRIX | apply_errmatrix) {
                        satband <- raster::raster(
                            (err_cube[,order_swir[band_swir], ]),
                            crs = "+proj=longlat +datum=WGS84")
                        satband <- raster::flip(satband, 1)
                        raster::projection(satband) <- NA
                    }
                    if (apply_errmatrix) {
                        band[satband > 0] <- NA
                    }
                }

            } else {
                if (proc_lev == "2D") {
                    message("Importing Band: ", band_swir,
                            " (", wl_swir[band_swir], ") of: of: 173")
                    outcrs <- paste0(
                        "+proj=utm +zone=", geo$proj_code,
                        ifelse(substring(
                            geo$proj_epsg, 3, 3) == 7, " +south", ""),
                        " +datum=WGS84 +units=m +no_defs")
                    band <- raster::raster(
                        (swir_cube[,order_swir[band_swir], ]),
                        crs = outcrs)
                    band <- raster::t(band)
                    band <- pr_setext_L2D(geo, band)
                    if (ERR_MATRIX | apply_errmatrix) {
                        satband <- raster::raster(
                            (err_cube[,order_swir[band_swir], ]),
                            crs = outcrs)
                        # satband <- raster::flip(satband, 1)
                        satband <- raster::t(satband)
                        satband <- pr_setext_L2D(geo, satband)
                    }
                    if (apply_errmatrix) {
                        band[satband > 0] <- NA
                    }
                }
            }
            # Add band to stack ----
            if (ind_band == 1) {
                rast_swir <- band
                if (ERR_MATRIX) rast_err <- satband
            } else {
                rast_swir <- raster::stack(rast_swir, band)
                if (ERR_MATRIX) rast_err <- raster::stack(rast_err, satband)
            }
            ind_band <- ind_band + 1
        } else {
            message("Band: ", band_swir, " not present")
        }
    }

    # Write the cube ----
    if (is.null(selbands_swir)) {
        orbands <- seqbands[wl_swir != 0]
        names(rast_swir) <- paste0("b", orbands)
        wl_sub   <- wl_swir[wl_swir != 0]
        fwhm_sub <- fwhm_swir[wl_swir != 0]
    } else {
        orbands <- seqbands
        names(rast_swir) <- paste0("b", orbands)
        wl_sub   <- wl_swir[seqbands]
        fwhm_sub <- fwhm_swir[seqbands]
    }
    rm(swir_cube)
    rm(band)
    gc()
    message("- Writing SWIR raster -")

    pr_rastwrite_lines(rast_swir,
                    out_file_swir,
                    out_format,
                    proc_lev,
                    scale_min = swir_min,
                    scale_max = swir_max)


    if (ERR_MATRIX) {
        message("- Writing ERR raster -")
        out_file_swir_err <- gsub("SWIR", "SWIR_ERR", out_file_swir)
        pr_rastwrite_lines(rast_err,
                        out_file_swir_err,
                        out_format,
                        "ERR",
                        scale_min = NULL,
                        scale_max = NULL)
        rm(rast_err)
    }

    if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file_swir), ".hdr")
        cat("band names = {", paste(names(rast_swir),collapse=","), "}", "\n",
            file=out_hdr, append=TRUE)
        write(c("wavelength = {",
                paste(round(wl_sub, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write(c("fwhm = {",
                paste(round(fwhm_sub, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write("wavelength units = Nanometers")
        write("sensor type = PRISMA")
    }

    rm(rast_swir)
    out_file_txt <- paste0(tools::file_path_sans_ext(out_file_swir), ".wvl")
    utils::write.table(data.frame(band = seq_along(wl_sub),
                                  orband = orbands,
                                  wl   = wl_sub,
                                  fwhm = fwhm_sub,
                                  stringsAsFactors = FALSE),
                       file = out_file_txt, row.names = FALSE)

}
