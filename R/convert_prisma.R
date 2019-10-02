#' @title convert_prisma
#' @description Access a PRISMA L1 HDF5 file and convert it to ENVI or GeoTiff
#'  format
#' @param in_file `character` full path of input HDF5 file
#' @param out_file `character` full path of output  file
#' @param out_format `character`` ["TIF" | "ENVI"], Output format, Default: 'tif'
#' @param source `character` ["HC0" | "HRC"], Considered Data Cube Default: 'HRC'
#' @param VNIR `logical` if TRUE, create the VNIR image, Default: TRUE
#' @param SWIR `logical` if TRUE, create the SWIR image, Default: TRUE
#' @param FULL `logical` if TRUE, create a single multispectral image from
#'  VNIR and SWIR, Default: FALSE
#' @param join_priority `character` ["VNIR" | "SWIR"], spectrometer to consider in
#'  the when join_spectra = TRUE, Default: SWIR - ignored if join_spectra is FALSE
#' @param PAN   `logical` if TRUE, also save the PAN data, default: TRUE
#' @param CLOUD `logical` if TRUE, also save the cloud mask data, default: TRUE
#' @param LC `logical` if TRUE, also save the land cover data, default: TRUE
#' @param overwrite `logical` if TRUE, existing files are overwritten, default: FALSE
#' @return The function is called for its side effects
#' @examples
#' \dontrun{
#'  in_file  <- "/home/lb/tmp/test/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
#'  out_file <- "/home/lb/tmp/test/test_1"
#'  out_format <- "ENVI"
#'
#'
#'  # Save a full image, prioritizing the SWIR spectrometer and save in EVI format
#'  convert_prisma(in_file       = in_file,
#'                 out_file      = out_file,
#'                 out_format    = out_format,
#'                 FULL          = TRUE,
#'                 join_priority = "SWIR",
#'                 LC            = TRUE,
#'                 CLOUD         = TRUE
#'                 )
#' }
#' @rdname convert_prisma
#' @export
#' @importFrom h5 h5file h5attr
#' @importFrom tools file_path_sans_ext
#' @importFrom raster stack raster t flip extent setExtent
#' @importFrom utils write.table
#'
convert_prisma <- function(in_file,
                           out_file,
                           out_format    = "ENVI",
                           VNIR          = TRUE,
                           SWIR          = TRUE,
                           FULL          = FALSE,
                           source        = "HRC",
                           join_priority = "SWIR",
                           PAN           = TRUE,
                           CLOUD         = FALSE,
                           LC            = FALSE,
                           overwrite     = FALSE) {

  # Open the file ----
  f <- try(h5::h5file(in_file))
  if (inherits(f, "try-error")) {
    stop("in_file does not appear to be a valid hdf5 dataset")
  }

  if (!dir.exists(dirname(out_file))) {
    if (dir.exists(dirname(dirname(out_file)))) {
      dir.create(dirname(out_file))
    } else {
      stop("Folder:", dirname(dirname(out_file)), " does not exist. Please create it beforehand!")
    }
  }

  # Get wavelengths and fwhms ----
  wl_vnir    <- h5::h5attr(f, "List_Cw_Vnir")
  order_vnir <- order(wl_vnir)
  wl_vnir <- wl_vnir[order_vnir]

  wl_swir    <- h5::h5attr(f, "List_Cw_Swir")
  order_swir <- order(wl_swir)
  wl_swir <- wl_swir[order_swir]

  wls <- c(wl_vnir, wl_swir)

  fwhm_vnir <- h5::h5attr(f, "List_Fwhm_Vnir")
  fwhm_vnir <- fwhm_vnir[order_vnir]

  fwhm_swir <- h5::h5attr(f, "List_Fwhm_Swir")
  fwhm_swir <- fwhm_swir[order_swir]

  fwhms <- c(fwhm_vnir, fwhm_swir)

  # get geolocation info ----
  lat <- t(f[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_SWIR")][])
  lon <- t(f[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_SWIR")][])

  # gett additional metadata
  sunzen  <- h5::h5attr(f, "Sun_zenith_angle")
  sunaz   <- h5::h5attr(f, "Sun_azimuth_angle")
  acqtime <- h5::h5attr(f, "Product_StartTime")

  out_file_angles <- paste0(tools::file_path_sans_ext(out_file), "_ANGLES.txt")
  utils::write.table(data.frame(date = acqtime,
                                sunzen   = sunzen,
                                sunaz = sunaz, stringsAsFactors = FALSE),
                     file = out_file_angles, row.names = FALSE)


  # get VNIR data cube and convert to raster ----

  out_file_vnir <- paste0(tools::file_path_sans_ext(out_file), "_VNIR")
  out_file_vnir <- ifelse(out_format == "TIF",
                          paste0(out_file_vnir, ".tif"),
                          paste0(out_file_vnir, ".envi"))

  if (VNIR) {
    message("- Importing VNIR Cube -")
    if (file.exists(out_file_vnir) & !overwrite) {
      message("VNIR file already exists - use overwrite = TRUE or change output file name to reprocess")
      rast_vnir <- raster::stack(out_file_vnir)
    } else {

      vnir_cube <- f[paste0("HDFEOS/SWATHS/PRS_L1_", source, "/Data Fields/VNIR_Cube")][]

      # vnir_cube <- aperm(vnir_cube, perm = c(1,3,2))[, , order_vnir]
      for (band_vnir in 1:66) {
        band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]), crs = "+proj=longlat +datum=WGS84")
        band <- raster::t(raster::flip(band, 2))
        ex <- matrix(c(min(lon), max(lon),  min(lat), max(lat)), nrow = 2, ncol = 2, byrow = T)
        ex <- raster::extent(ex)
        band <- raster::setExtent(band, ex, keepres=F)
        if (band_vnir == 1) {
          rast_vnir <- band
        } else {
          rast_vnir <- raster::stack(rast_vnir, band)
        }
      }
      rm(vnir_cube)
      rm(band)
      gc()

      if (VNIR) {
        message("- Writing VNIR raster -")
        rastwrite_lines(rast_vnir, out_file_vnir, out_format)
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
    }
  }
  # get SWIR data cube and convert to raster ----
  #
  out_file_swir <- paste0(tools::file_path_sans_ext(out_file), "_SWIR")
  out_file_swir <- ifelse(out_format == "TIF",
                          paste0(out_file_swir, ".tif"),
                          paste0(out_file_swir, ".envi"))
  if (SWIR) {
    if (file.exists(out_file_swir) & !overwrite) {
      message("SWIR file already exists - use overwrite = TRUE or change output file name to reprocess")
      rast_swir <- raster::stack(out_file_swir)
    } else {

      message("- Importing SWIR Cube - ")
      swir_cube <- f[paste0("HDFEOS/SWATHS/PRS_L1_", source, "/Data Fields/SWIR_Cube")][]
      # swir_cube <- aperm(swir_cube, perm = c(1,3,2))[, , order_swir]
      for (band_swir in 1:173) {
        band <- raster::raster((swir_cube[,order_swir[band_swir], ]), crs = "+proj=longlat +datum=WGS84")
        band <- raster::t(raster::flip(band, 2))
        ex   <- matrix(c(min(lon), max(lon),  min(lat), max(lat)), nrow = 2, ncol = 2, byrow = T)
        ex   <- raster::extent(ex)
        band <- raster::setExtent(band, ex, keepres=F)
        if (band_swir == 1) {
          rast_swir <- band
        } else {
          rast_swir <- raster::stack(rast_swir, band)
        }
      }
      rm(swir_cube)
      rm(band)
      gc()

      if (SWIR) {
        message("- Writing SWIR raster -")
        rastwrite_lines(rast_swir, out_file_swir, out_format)
        if (out_format == "ENVI") {
          out_hdr <- paste0(tools::file_path_sans_ext(out_file_swir), ".hdr")

          # writeLines(c("wavelength units = DOY"), fileConn_meta_hdr)
          # Wavelengths == DOY from 01/01/2000
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
    }
  }

  # create FULL data cube and convert to raster ----
  out_file_full <- paste0(tools::file_path_sans_ext(out_file), "_FULL")
  out_file_full <- ifelse(out_format == "TIF",
                          paste0(out_file_full, ".tif"),
                          paste0(out_file_full, ".envi"))
  if (FULL) {

    if (file.exists(out_file_full) & !overwrite) {
      message("FULL file already exists - use overwrite = TRUE or change output file name to reprocess")
    } else {
      message("- Creating FULL raster -")
      # Save hyperspectral cube
      if (join_priority == "SWIR") {
        rast_tot <- raster::stack(rast_vnir[[1:58]], rast_swir)
        wl_tot   <- c(wl_vnir[1:58], wl_swir)
        fwhm_tot <- c(fwhm_vnir[1:58], fwhm_swir)
      } else {
        rast_tot <- raster::stack(rast_vnir, rast_swir[10:173])
        wl_tot   <- c(wl_vnir, wl_swir[10:173])
        fwhm_tot <- c(fwhm_vnir, fwhm_swir[10:173])
      }

      rm(rast_vnir)
      rm(rast_swir)
      gc()

      message("- Writing FULL raster -")
      rastwrite_lines(rast_tot, out_file_full, out_format)

      if (out_format == "ENVI") {
        out_hdr <- paste0(tools::file_path_sans_ext(out_file_full), ".hdr")

        # writeLines(c("wavelength units = DOY"), fileConn_meta_hdr)
        # Wavelengths == DOY from 01/01/2000
        write(c("wavelength = {",
                paste(round(wl_tot, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
        write(c("fwhm = {",
                paste(round(fwhm_tot, digits = 4), collapse = ","), "}"),
              out_hdr, append = TRUE)
      }
      out_file_txt <- paste0(tools::file_path_sans_ext(out_file_full), "_meta.txt")
      utils::write.table(data.frame(band = 1:length(wl_tot),
                                    wl   = wl_tot,
                                    fwhm = fwhm_tot, stringsAsFactors = FALSE),
                         file = out_file_txt, row.names = FALSE)

      rm(rast_tot)
      gc()
    }

  }

  # Save PAN if requested ----
  out_file_pan <- paste0(tools::file_path_sans_ext(out_file), "_PAN")
  out_file_pan <- ifelse(out_format == "TIF",
                         paste0(out_file_pan, ".tif"),
                         paste0(out_file_pan, ".envi"))
  if (file.exists(out_file_pan) & !overwrite) {
    message("PAN file already exists - use overwrite = TRUE or change output file name to reprocess")
  } else {

    if (PAN) {

      message(" - Accessing PAN raster - ")

      pan_cube <- f[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source), "/Data Fields/Cube")][]
      pan_lat <- t(f[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                            "/Geolocation Fields/Latitude")][])
      pan_lon <- t(f[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                            "/Geolocation Fields/Longitude")][])
      rast_pan <- raster::raster(pan_cube, crs = "+proj=longlat +datum=WGS84")
      rast_pan <- raster::t(raster::flip(rast_pan, 2))
      rm(pan_cube)
      gc()

      ex   <- matrix(c(min(pan_lon), max(pan_lon),
                       min(pan_lat), max(pan_lat)),
                     nrow = 2, ncol = 2, byrow = T)
      ex   <- raster::extent(ex)
      rast_pan <- raster::setExtent(rast_pan, ex, keepres=F)

      message("- Writing PAN raster -")

      rastwrite_lines(rast_pan, out_file_pan, out_format)
      rm(rast_pan)
      rm(pan_lon)
      rm(pan_lat)
      gc()
    }
  }

  # Save CLD if requested ----
  out_file_cld <- paste0(tools::file_path_sans_ext(out_file), "_CLD")
  out_file_cld <- ifelse(out_format == "TIF",
                         paste0(out_file_cld, ".tif"),
                         paste0(out_file_cld, ".envi"))

  if (file.exists(out_file_cld) & !overwrite) {
    message("CLD file already exists - use overwrite = TRUE or change output file name to reprocess")
  } else {
    if (CLOUD) {

      message(" - Accessing CLOUD raster - ")
      cld_cube <- f["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/Cloud_Mask"][]
      rast_cld <- raster::raster(cld_cube, crs = "+proj=longlat +datum=WGS84")
      rast_cld <- raster::t(raster::flip(rast_cld, 2))
      rm(cld_cube)
      gc()

      ex   <- matrix(c(min(lon), max(lon),
                       min(lat), max(lat)),
                     nrow = 2, ncol = 2, byrow = T)
      ex   <- raster::extent(ex)
      rast_cld <- raster::setExtent(rast_cld, ex, keepres=F)

      message("- Writing CLD raster -")
      rastwrite_lines(rast_cld, out_file_cld, out_format)

    }
  }

  # Save LC if requested ----
  out_file_lc <- paste0(tools::file_path_sans_ext(out_file), "_LC")
  out_file_lc <- ifelse(out_format == "TIF",
                        paste0(out_file_lc, ".tif"),
                        paste0(out_file_lc, ".envi"))

  if (LC) {
    if (file.exists(out_file_lc) & !overwrite) {
      message("LC file already exists - use overwrite = TRUE or change output file name to reprocess")
    } else {

      message(" - Accessing LAND COVER raster - ")
      lc_cube <- f["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/LandCover_Mask"][]
      rast_lc <- raster::raster(lc_cube, crs = "+proj=longlat +datum=WGS84")
      rast_lc <- raster::t(raster::flip(rast_lc, 2))
      rm(lc_cube)
      gc()

      ex   <- matrix(c(min(lon), max(lon),
                       min(lat), max(lat)),
                     nrow = 2, ncol = 2, byrow = T)
      ex   <- raster::extent(ex)
      rast_lc <- raster::setExtent(rast_lc, ex, keepres=F)

      message("- Writing LC raster -")

      rastwrite_lines(rast_lc, out_file_lc, out_format)

    }

  }

}
