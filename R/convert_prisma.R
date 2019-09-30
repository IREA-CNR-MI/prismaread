#' @title convert_prisma
#' @description Access a PRISMA L1 HDF5 file and convert it to ENVI or GeoTiff
#'  format
#' @param in_file `character` full path of input HDF5 file
#' @param out_file `character` full path of output  file
#' @param out_format `character`` ["tif" | "ENVI"], Output format, Default: 'tif'
#' @param source `character` ["HC0" | "HRC"], Considered Data Cube Default: 'HRC'
#' @param join_spectra `logical` if TRUE, create a single multispectral image from
#'  VNIR and SWIR, otherwise, save two separate images, Default: FALSE
#' @param join_priority `character` ["VNIR" | "SWIR"], spectrometer to consider in
#'  the when join_spectra = TRUE, Default: SWIR - ignored if join_spectra is FALSE
#' @param pan   `logical` if TRUE, also save the PAN data, default: TRUE
#' @param cloud `logical` if TRUE, also save the cloud mask data, default: TRUE
#' @param landcov `logical` if TRUE, also save the land cover data, default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  in_file  <- "/home/lb/projects/ASI-PRISCAV/3_IMAGES/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
#'  out_file <- "/home/lb/projects/test/test_1_all.envi"
#'
#'  # Save a full image, prioritizing the SWIR spectrometer and save in EVI format
#'  convert_prisma(in_file       = in_file,
#'                 out_file      = out_file,
#'                 format        = "ENVI",
#'                 join_spectra  = TRUE,
#'                 join_priority = "SWIR"
#'                 )
#'  }
#' }
#' @seealso
#'  \code{\link[h5]{c("H5File", "H5File")}},\code{\link[h5]{H5Location-Attribute}}
#'  \code{\link[raster]{c("raster", "Raster")}},\code{\link[raster]{transpose}},\code{\link[raster]{c("flip", "flip")}},\code{\link[raster]{c("Extent-class", "extent", "extent")}},\code{\link[raster]{setExtent}},\code{\link[raster]{stack}},\code{\link[raster]{brick}},\code{\link[raster]{blockSize}},\code{\link[raster]{writeValues}},\code{\link[raster]{getValues}}
#'  \code{\link[tools]{fileutils}}
#'  \code{\link[mapview]{c("mapView", "mapView")}},\code{\link[mapview]{viewRGB}}
#' @rdname convert_prisma
#' @export
#' @importFrom h5 h5file h5attr
#' @importFrom raster raster t flip extent setExtent stack brick blockSize writeStart getValues writeValues writeStop
#' @importFrom tools file_path_sans_ext
#' @importFrom mapview mapview viewRGB
#'
convert_prisma <- function(in_file,
                           out_file,
                           out_format = "tif",
                           source     = "HRC",
                           join_spectra = FALSE) {

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


  # get VNIR data cube and convert to raster ----

  message("- Importing VNIR Cube -")

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

  if (join_spectra == FALSE) {

    out_file <- paste(tools::file_path_sans_ext(out_file), "_SWIR")
    out_file <- ifelse(out_format == "tif",
                       paste0(out_file, ".tif"),
                       paste0(out_file, ".envi"))
    rastwrite_lines(rast_vnir, out_file, out_format)
    if (out_format == "ENVI") {
      out_hdr <- paste0(tools::file_path_sans_ext(out_file), ".hdr")

      # writeLines(c("wavelength units = DOY"), fileConn_meta_hdr)
      # Wavelengths == DOY from 01/01/2000
      write(c("wavelength = {",
              paste(round(wls_vnir, digits = 4), collapse = ","), "}"),
            out_hdr, append = TRUE)
    }

  }

  # get SWIR data cube and convert to raster ----
  #

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

  if (join_spectra == FALSE) {

    out_file <- paste(tools::file_path_sans_ext(out_file), "_SWIR")
    out_file <- ifelse(out_format == "tif",
                       paste0(out_file, ".tif"),
                       paste0(out_file, ".envi"))
    rastwrite_lines(rast_swir, out_file, out_format)
    if (out_format == "ENVI") {
      out_hdr <- paste0(tools::file_path_sans_ext(out_file), ".hdr")

      # writeLines(c("wavelength units = DOY"), fileConn_meta_hdr)
      # Wavelengths == DOY from 01/01/2000
      write(c("wavelength = {",
              paste(round(wls_vnir, digits = 4), collapse = ","), "}"),
            out_hdr, append = TRUE)
    }

  }


  if (join_spectra == TRUE) {

    # Save hyperspectral cube
    rast_tot <- raster::stack(rast_vnir, rast_swir)
    rm(rast_vnir)
    rm(rast_swir)
    gc()
    message("- Saving SWIR Cube - ")

    if (format == "tif") {
      out <- raster::brick(rast_tot, values = FALSE)
      bs <-  raster::blockSize(out)
      out <- raster::writeStart(out, filename = out_file, overwrite = TRUE)

      for (i in 1:bs$n) {
        print(i)
        v <- raster::getValues(rast_tot, row=bs$row[i], nrows=bs$nrows[i] )
        out <- raster::writeValues(out, v, bs$row[i])
      }
      out <- raster::writeStop(out)
    } else {
      out <- raster::brick(rast_tot, values = FALSE)
      bs <-  raster::blockSize(out)
      out <- raster::writeStart(out, filename = out_file, overwrite = TRUE,
                                format = "ENVI")

      for (i in 1:bs$n) {
        print(i)
        v <- raster::getValues(rast_tot, row=bs$row[i], nrows=bs$nrows[i] )
        out <- raster::writeValues(out, v, bs$row[i])
      }

      out <- raster::writeStop(out)
      out_hdr <- paste0(tools::file_path_sans_ext(out_file), ".hdr")

      # writeLines(c("wavelength units = DOY"), fileConn_meta_hdr)
      # Wavelengths == DOY from 01/01/2000
      write(c("wavelength = {",
              paste(round(wls, digits = 4), collapse = ","), "}"),
            out_hdr, append = TRUE)
    }

  }

  # stitch vnir and swir and save

  # return(out)

  gc()

  pan_Cube <- f["/HDFEOS/SWATHS/PRS_L1_PRC/Data Fields/Cube"][]
  aa <- setExtent(aa, ex, keepres=F)






  # mat_tot <- c(mat_vnir, mat_swir)






  rm(rast_vnir)
  rm(rast_swir)
  rm(mat_vnir)
  rm(mat_swir)


  crs(rast_tot) <- "+proj=longlat +datum=WGS84"
  mapview::mapview(rast_tot[[1]])

  mapview::viewRGB(rast_tot, 40,30, 20)

  spect <- rast_tot[100,100][1,]
  dfp <- data.frame(wl = wls, Rad =spect)
  ggplot(dfp, aes(x = wl, y = Rad)) + geom_line() + theme_light()


  # writeRaster("")
}
