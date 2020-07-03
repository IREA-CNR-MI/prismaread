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
                               fill_gaps,
                               fix_geo,
                               wl_swir,
                               order_swir,
                               fwhm_swir,
                               apply_errmatrix,
                               ERR_MATRIX,
                               selbands_swir = NULL,
                               in_L2_file = NULL){

    # Get geo info ----
    geo <- prisma_get_geoloc(f, proc_lev, source, wvl = "SWIR", in_L2_file)

    # Get the datacube and required attributes frim hdr ----
    if (proc_lev == "1") {
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source,
                               "/Data Fields/SWIR_Cube")]][,,]
        swir_scale  <- hdf5r::h5attr(f, "ScaleFactor_Swir")
        swir_offset <- hdf5r::h5attr(f, "Offset_Swir")
    } else {
        swir_max  <- hdf5r::h5attr(f, "L2ScaleSwirMax")
        swir_min  <- hdf5r::h5attr(f, "L2ScaleSwirMin")
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L", proc_lev,"_",
                               source, "/Data Fields/SWIR_Cube")]][,,]
    }

    # Get the different bands in order of wvl, and convert to `raster` bands ----
    # Also georeference if needed
    ind_band <- 1
    if (is.null(selbands_swir)) {
        seqbands <- 1:173
    } else {
        seqbands <- unlist(lapply(selbands_swir, FUN = function(x) which.min(abs(wl_swir - x))))
    }

    for (band_swir in seqbands) {

        if (wl_swir[band_swir] != 0) {
            if (proc_lev %in% c("1", "2B", "2C")) {
                if (base_georef) {
                    message("Importing Band: ", band_swir, " of: 173 and applying bowtie georeferencing")
                    band <- raster::raster((swir_cube[,order_swir[band_swir], ]),
                                           crs = "+proj=longlat +datum=WGS84")
                    lat  <- geo$lat
                    lon  <- geo$lon
                    if (proc_lev == "1") {
                        band <- (band / swir_scale) - swir_offset
                    }
                    band <- prisma_basegeo(band, lon, lat, fill_gaps)
                } else {
                    message("Importing Band: ", band_swir, " of: 173")
                    band <- raster::raster(swir_cube[,order_swir[band_swir], ],
                                           crs = "+proj=longlat +datum=WGS84")
                    if (proc_lev == "1") {
                        band <- (band / swir_scale) - swir_offset
                    }
                    band <- raster::flip(band, 1)
                }

            } else {
                if (proc_lev == "2D") {
                    message("Importing Band: ", band_swir, " of: 173")
                    band <- raster::raster((swir_cube[,order_swir[band_swir], ]),
                                           crs = paste0("+proj=utm +zone=", geo$proj_code,
                                                        ifelse(substring(geo$proj_epsg, 3, 3) == 7, " +south", ""),
                                                        " +datum=WGS84 +units=m +no_defs"))
                    band <- raster::t(band)
                    ex <- matrix(c(geo$xmin - 15, geo$xmin - 15 + dim(band)[2]*30,
                                   geo$ymin - 15, geo$ymin - 15 + dim(band)[1]*30),
                                 nrow = 2, ncol = 2, byrow = T)
                    ex <- raster::extent(ex)
                    band <- raster::setExtent(band, ex, keepres = FALSE)
                }
            }
            # Add band to stack ----
            if (ind_band == 1) {
                rast_swir <- band
            } else {
                rast_swir <- raster::stack(rast_swir, band)
            }
            ind_band <- ind_band + 1
        } else {
            message("Band: ", band_swir, " not present")
        }
    }

    # Write the cube ----
    if (is.null(selbands_swir)) {
        # names(rast_vnir) <- paste0("b", seqbands, "_", round(wl_vnir, digits = 3))
        orbands <- seqbands[wl_swir != 0]
        names(rast_swir) <- paste0("b", orbands)
        wl_sub   <- wl_swir[wl_swir != 0]
        fwhm_sub <- fwhm_swir[wl_swir != 0]
    } else {
        # names(rast_vnir) <- paste0("b", seqbands, "_", round(wl_vnir[seqbands], digits = 3))
        orbands <- seqbands
        names(rast_swir) <- paste0("b", orbands)
        wl_sub   <- wl_swir[seqbands]
        fwhm_sub <- fwhm_swir[seqbands]
    }
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
        cat("band names = {", paste(names(rast_swir),collapse=","), "}", "\n",
            file=raster::extension(out_file_swir, "hdr"), append=TRUE)
        out_hdr <- paste0(tools::file_path_sans_ext(out_file_swir), ".hdr")
        write(c("wavelength = {",
                paste(round(wl_sub, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write(c("fwhm = {",
                paste(round(fwhm_sub, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write("wavelength units = Nanometers")
        write("sensor type = PRISMA")
    }

    out_file_txt <- paste0(tools::file_path_sans_ext(out_file_swir), ".wvl")
    utils::write.table(data.frame(band = 1:length(wl_sub),
                                  orband = orbands,
                                  wl   = wl_sub,
                                  fwhm = fwhm_sub,
                                  stringsAsFactors = FALSE),
                       file = out_file_txt, row.names = FALSE)

}
