#' @title pr_create_vnir
#' @description helper function used to process and save the VNIR data cube
#' @param f input data he5 from caller
#' @param proc_lev `character` Processing level (e.g., "1", "2B") - passed by
#'  caller
#' @param out_file_vnir output file name for VNIR
#' @param wl_vnir passed by caller - array of PRISMA VNIR wavelengths
#' @param order_vnir passed by caller - ordering of array of PRISMA VNIR
#'   wavelengths
#' @param fwhm_vnir passed by caller - array of PRISMA VNIR fwhms
#' @inheritParams pr_convert
#' @return the function is called for its side effects
#' @importFrom hdf5r h5attr
#' @importFrom raster raster extent flip t setExtent stack
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#'
pr_create_vnir <- function(f,
                               proc_lev,
                               source,
                               out_file_vnir,
                               out_format,
                               base_georef,
                               fill_gaps,
                               wl_vnir,
                               order_vnir,
                               fwhm_vnir,
                               apply_errmatrix,
                               ERR_MATRIX,
                               selbands_vnir = NULL,
                               in_L2_file = NULL){

    # Get geo info ----

    geo <- pr_get_geoloc(f, proc_lev, source, wvl = "VNIR", in_L2_file)

    # Get the datacube and required attributes from hdr ----
    if(proc_lev == 1) {
        vnir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source,
                               "/Data Fields/VNIR_Cube")]][,,]
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
            err_cube <- f[[paste0(
                "HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                source, "/Data Fields/VNIR_PIXEL_L2_ERR_MATRIX")]][,,]
        }
    }

    # Get the different bands in order of wvl, and convert to `raster` bands----
    # Also georeference if needed
    ind_vnir <- 1
    if (is.null(selbands_vnir)) {
        seqbands <- 1:66
    } else {
        seqbands <- unlist(lapply(
            selbands_vnir, FUN = function(x) which.min(abs(wl_vnir - x))))
    }

    for (band_vnir in seqbands) {

        if (wl_vnir[band_vnir] != 0) { #skip 0-wavelength bands

            if(proc_lev %in% c("1", "2B", "2C")) {
                # on L1, 2B or 2C, apply bowtie georeferencing if requested ----
                band <- raster::raster((vnir_cube[,order_vnir[band_vnir],]),
                                       crs = "+proj=longlat +datum=WGS84")
                if (base_georef) {
                    message("Importing Band: ", band_vnir,
                            " (",wl_vnir[band_vnir], ") of: 66 and
                            applying bowtie georeferencing")
                    lat  <- geo$lat
                    lon  <- geo$lon
                    if (proc_lev == "1") {
                        band <- (band / vnir_scale) - vnir_offset
                    }
                    band <- pr_basegeo(band, lon, lat, fill_gaps)
                    if (apply_errmatrix | ERR_MATRIX) {
                        satband <- raster::raster(
                            (err_cube[,order_vnir[band_vnir], ]),
                            crs = "+proj=longlat +datum=WGS84")
                        satband <- pr_basegeo(satband, lon, lat, fill_gaps)
                    }
                    if (apply_errmatrix) {
                        band[satband > 0] <- NA
                    }

                } else {
                    message("Importing Band: ", band_vnir,
                            " (",wl_vnir[band_vnir], ") of: 66")
                    if (proc_lev == "1") {
                        band <- (band / vnir_scale) - vnir_offset
                    }
                    # flip the band to get it north/south
                    band <- raster::flip(band, 1)
                    raster::projection(band) <- NA
                    if (apply_errmatrix | ERR_MATRIX) {
                        satband <- raster::raster(
                            (err_cube[,order_vnir[band_vnir], ]),
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
                    # on L2, retreive and apply georeferencing ----
                    message("Importing Band: ", band_vnir,
                            " (",wl_vnir[band_vnir], ") of: 66")
                    outcrs <- paste0(
                        "+proj=utm +zone=", geo$proj_code,
                        ifelse(substring(
                            geo$proj_epsg, 3, 3) == 7, " +south", ""),
                        " +datum=WGS84 +units=m +no_defs")
                    band <- raster::raster(
                        (vnir_cube[,order_vnir[band_vnir], ]),
                        crs = outcrs)
                    # traspose the band to get it correctly oriented and set
                    # extent
                    band <- raster::t(band)
                    band <- pr_setext_L2D(geo, band)
                    if (apply_errmatrix | ERR_MATRIX) {
                        satband <- raster::raster(
                            err_cube[,order_vnir[band_vnir], ],
                            crs = outcrs)
                        satband <- raster::t(satband)
                        satband <- pr_setext_L2D(geo, satband)
                    }
                    if (apply_errmatrix) {
                        band[satband > 0] <- NA
                    }

                }
            }

            # grow the raster stack by adding current band ----
            if (ind_vnir == 1) {
                rast_vnir <- band
                if (ERR_MATRIX) rast_err <- satband
            } else {
                rast_vnir <- raster::stack(rast_vnir, band)
                if (ERR_MATRIX) rast_err <- raster::stack(rast_err, satband)
            }
            ind_vnir <- ind_vnir + 1

        } else {
            message("Band: ", band_vnir, " not present")
        }
    }

    # Write the cube ----

    # create arrays of wavelengths to be used for creation of envi header and
    # bandnames
    if (is.null(selbands_vnir)) {
        orbands <- seqbands[wl_vnir != 0]
        names(rast_vnir) <- paste0("b", orbands)
        wl_sub   <- wl_vnir[wl_vnir != 0]
        fwhm_sub <- fwhm_vnir[wl_vnir != 0]
    } else {
        orbands <- seqbands
        names(rast_vnir) <- paste0("b", orbands)
        wl_sub   <- wl_vnir[seqbands]
        fwhm_sub <- fwhm_vnir[seqbands]
    }

    rm(vnir_cube)
    rm(band)
    gc()

    message("- Writing VNIR raster -")

    pr_rastwrite_lines(rast_vnir,
                    out_file_vnir,
                    out_format,
                    proc_lev,
                    scale_min = vnir_min,
                    scale_max = vnir_max)


    if (ERR_MATRIX) {
        message("- Writing ERR raster -")
        out_file_vnir_err <- gsub("VNIR", "VNIR_ERR", out_file_vnir)
        pr_rastwrite_lines(rast_err,
                        out_file_vnir_err,
                        out_format,
                        "ERR",
                        scale_min = NULL,
                        scale_max = NULL)
        rm(rast_err)
    }


    if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file_vnir), ".hdr")
        cat("band names = {", paste(names(rast_vnir),collapse=","), "}", "\n",
            file=raster::extension(out_file_vnir, "hdr"), append=TRUE)
        write(c("wavelength = {",
                paste(round(wl_sub, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write(c("fwhm = {",
                paste(round(fwhm_sub, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write("wavelength units = Nanometers")
        write("sensor type = PRISMA")

    }
    rm(rast_vnir)
    # write textfile of wavelengths ----
    out_file_txt <- paste0(tools::file_path_sans_ext(out_file_vnir), ".wvl")
    utils::write.table(data.frame(band = seq_along(wl_sub),
                                  orband = orbands,
                                  wl   = wl_sub,
                                  fwhm = fwhm_sub,
                                  stringsAsFactors = FALSE),
                       file = out_file_txt, row.names = FALSE)
}
