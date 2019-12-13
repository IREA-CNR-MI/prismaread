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
#'  the when join_spectra = TRUE, Default: SWIR - ignored if join_spectra is FALSE.
#'  Default: "VNIR"
#' @param ATCOR  logical` if TRUE, create the text files required to run ATCOR, Default: FALSE;
#' @param ATCOR_wls  `character`, or `numeric` If equal to `Nominal`, the only ATCOR wvl file
#'  created is the one containing Nominal wavelengths. If a numeric vector is provided,
#'  then one different wvl file is created for each selected COLUMN selected (e.g., if
#'  providing `ATCOR_wls = c(200, 800)`, then the wavelengths and FWHMs related to
#'  columns 200 and 800 are saved.)
#' @param PAN   `logical` if TRUE, also save the PAN data, default: TRUE
#' @param CLOUD `logical` if TRUE, also save the cloud mask data, default: TRUE
#' @param GLINT `logical` if TRUE, also save the GLINT mask data, default: TRUE
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
#'  # Save a full image, prioritizing the VNIR spectrometer and save in ENVI format
#'  convert_prisma(in_file       = in_file,
#'                 out_file      = out_file,
#'                 out_format    = out_format,
#'                 FULL          = TRUE,
#'                 join_priority = "VNIR",
#'                 LC            = TRUE,
#'                 CLOUD         = TRUE,
#'                 overwrite     = TRUE
#'                 )
#'
#'  # Save a full image, prioritizing the SWIR spectrometer and save in ENVI format,
#'  # Also create ATCOR files, with Nominal wavelengths, and those for columns
#'  # 200 and 800 of the cube
#'  convert_prisma(in_file       = in_file,
#'                 out_file      = out_file,
#'                 out_format    = out_format,
#'                 FULL          = TRUE,
#'                 join_priority = "VNIR",
#'                 ATCOR         = TRUE,
#'                 ATCOR_wls     = c(200,800),
#'                 LC            = TRUE,
#'                 CLOUD         = TRUE,
#'                 overwrite     = TRUE
#'                 )
#'
#'
#' }
#' @rdname convert_prisma
#' @export
#' @importFrom hdf5r H5File h5attr
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
                           source        = "HCO",
                           join_priority = "SWIR",
                           ATCOR         = TRUE,
                           ATCOR_wls    = "Nominal",
                           PAN           = TRUE,
                           CLOUD         = FALSE,
                           LC            = FALSE,
                           GLINT         = FALSE,
                           overwrite     = FALSE) {

  # browser()

  # Open the file ----
  f <- hdf5r::H5File$new(in_file, mode="r+")
  proc_lev <- hdf5r::h5attr(f, "Processing_Level")

  if (proc_lev != "1") {
    if (source %in% c("HRC", "HC0")) {
      message("Processing Level = 2 - Source modified to \"HCO\" by default")
      source = "HCO"
    }
  }


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
  wl_vnir    <- hdf5r::h5attr(f, "List_Cw_Vnir")
  order_vnir <- order(wl_vnir)
  wl_vnir <- wl_vnir[order_vnir]

  wl_swir    <- hdf5r::h5attr(f, "List_Cw_Swir")
  order_swir <- order(wl_swir)
  wl_swir <- wl_swir[order_swir]
  wls <- c(wl_vnir, wl_swir)

  fwhm_vnir <- hdf5r::h5attr(f, "List_Fwhm_Vnir")
  fwhm_vnir <- fwhm_vnir[order_vnir]

  fwhm_swir <- hdf5r::h5attr(f, "List_Fwhm_Swir")
  fwhm_swir <- fwhm_swir[order_swir]

  fwhms <- c(fwhm_vnir, fwhm_swir)

  # write ATCOR files if needed ----
  if (ATCOR == TRUE && proc_lev == "1") {

    ATCOR_fold <- file.path(dirname(out_file), "ATCOR")
    dir.create(ATCOR_fold, showWarnings = FALSE)
    out_file_wvl <- file.path(ATCOR_fold,
                              paste0(tools::file_path_sans_ext(basename(out_file)),
                                     "_atcor_wvl_nominal.wvl"))
    wl_tot_atcor   <- wls[wls != 0]
    fwhm_tot_atcor <- fwhms[fwhms != 0]
    out <- data.frame(`channel number`            = 1:length(wl_tot_atcor),
                      `channel center wavelength` = round(wl_tot_atcor/1000, digits = 6),
                      `bandwidth` = fwhm_tot_atcor, stringsAsFactors = FALSE)
    names(out) <- c("channel number", "channel center wavelength", "bandwidth")
    utils::write.table(out, file = out_file_wvl, row.names = FALSE, sep = "\t")

    out_file_cal <- file.path(ATCOR_fold,
                              paste0(tools::file_path_sans_ext(basename(out_file)),
                                     "_atcor_cal.cal"))
    out <- data.frame("wavelength" = round(wl_tot_atcor/1000, digits = 6),
                      `radiometric offset c0` = 0,
                      `radiometric offset c1` = 1, stringsAsFactors = FALSE)
    names(out) <- c("wavelength", "radiometric offset c0 ", "radiometric offset c0")
    utils::write.table(out, file = out_file_cal, row.names = FALSE, sep = "\t")

    out_file_dat <- file.path(ATCOR_fold,
                              paste0(tools::file_path_sans_ext(basename(out_file)),
                                     "atcor_dat.dat"))
    file.copy(system.file("extdata/atcor_dat.dat", package = "prismaread"),
              out_file_dat)

    # if specified, save additional wvl files corresponding to selected columns ----

    if (ATCOR_wls != "Nominal" && ATCOR == TRUE) {

      if (hdf5r::existsGroup(f, "//KDP_AUX/Cw_Vnir_Matrix")) {

        if (!is.numeric(ATCOR_wls)) {
          stop("ATCOR_wls should be wither \"Nominal\" or a vector containing the",
               "column numbers at which wavelengths should be retrieved")
        }

        vnir_start  <- hdf5r::h5attr(f, "Start_index_EO_VNIR")
        vnir_stop   <- hdf5r::h5attr(f, "Stop_index_EO_VNIR")
        vnir_wl_mat <- t(f[["//KDP_AUX/Cw_Vnir_Matrix"]][vnir_start:vnir_stop,])
        vnir_wl_mat <- vnir_wl_mat[,order_vnir]
        vnir_wl_mat <- vnir_wl_mat[, which(vnir_wl_mat[1,] != 0)]
        vnir_fwhm_mat <- t(f[["KDP_AUX/Fwhm_Vnir_Matrix"]][vnir_start:vnir_stop,])
        vnir_fwhm_mat <- vnir_fwhm_mat[,order_vnir]
        vnir_fwhm_mat <- vnir_fwhm_mat[, which(vnir_fwhm_mat[1,] != 0)]

        swir_start  <- hdf5r::h5attr(f, "Start_index_EO_SWIR")
        swir_stop   <- hdf5r::h5attr(f, "Stop_index_EO_SWIR")
        swir_wl_mat <- t(f[["//KDP_AUX/Cw_Swir_Matrix"]][swir_start:swir_stop,])
        swir_wl_mat <- swir_wl_mat[,order_swir]
        swir_wl_mat <- swir_wl_mat[, which(swir_wl_mat[1,] != 0)]
        swir_fwhm_mat <- t(f[["//KDP_AUX/Fwhm_Swir_Matrix"]][swir_start:swir_stop,])
        swir_fwhm_mat <- swir_fwhm_mat[,order_swir]
        swir_fwhm_mat <- swir_fwhm_mat[, which(swir_fwhm_mat[1,] != 0)]

        if(join_priority == "VNIR") {
          swir_wl_mat   <- swir_wl_mat[,which(swir_wl_mat[1,] > max(vnir_wl_mat[1,]))]
          swir_fwhm_mat <- swir_fwhm_mat[,which(swir_wl_mat[1,] > max(vnir_wl_mat[1,]))]
        } else {
          vnir_wl_mat   <- vnir_wl_mat[,which(vnir_wl_mat[1,] < min(swir_wl_mat[1,]))]
          vnir_fwhm_mat <- vnir_fwhm_mat[,which(vnir_wl_mat[1,] < min(swir_wl_mat[1,]))]
        }

        wl_mat_tot    <- cbind(vnir_wl_mat, swir_wl_mat )
        fwhm_mat_tot  <- cbind(vnir_fwhm_mat, swir_fwhm_mat )

        for (col in ATCOR_wls) {
          dir.create(file.path(ATCOR_fold, trimws(col)), showWarnings = FALSE)
          out_file_wvl <- file.path(ATCOR_fold, trimws(col),
                                    paste0(tools::file_path_sans_ext(basename(out_file)),
                                           paste0("_atcor_wvl_", trimws(col), ".wvl")))
          out <- data.frame(`channel number`            = 1:dim(wl_mat_tot)[2],
                            `channel center wavelength` = round(wl_mat_tot[col,]/1000, digits = 6),
                            `bandwidth` = round(fwhm_mat_tot[col,], digits = 6),
                            stringsAsFactors = FALSE)
          names(out) <- paste(c("channel number", "channel center wavelength", "bandwidth"), col)
          utils::write.table(out, file = out_file_wvl, row.names = FALSE, sep = "\t")

          out_file_cal <- file.path(ATCOR_fold, trimws(col),
                                    paste0(tools::file_path_sans_ext(basename(out_file)),
                                           "_atcor_cal.cal"))
          out <- data.frame("wavelength" = round(wl_mat_tot[col,]/1000, digits = 6),
                            `radiometric offset c0` = 0,
                            `radiometric offset c1` = 1, stringsAsFactors = FALSE)
          names(out) <- c("wavelength", "radiometric offset c0 ", "radiometric offset c0")
          utils::write.table(out, file = out_file_cal, row.names = FALSE, sep = "\t")

          out_file_dat <- file.path(ATCOR_fold, trimws(col),
                                    paste0(tools::file_path_sans_ext(basename(out_file)),
                                           "atcor_dat.dat"))
          file.copy(system.file("extdata/atcor_dat.dat", package = "prismaread"),
                    out_file_dat)
        }
      } else {
        message("CW matrix dataset not existing - creation of additional ATCOR files ignored")
      }
    }
  }

  # get geolocation info ----
  if (proc_lev == "1") {
    lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Latitude_SWIR")]][,])
    lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", source, "/Geolocation Fields/Longitude_SWIR")]][,])
  } else {
    proj_code <- hdf5r::h5attr(f, "Projection_Id")
    proj_name <- hdf5r::h5attr(f, "Projection_Name")
    xmin  <- hdf5r::h5attr(f, "Product_ULcorner_easting")
    xmax  <- hdf5r::h5attr(f, "Product_LRcorner_easting")
    ymin  <- hdf5r::h5attr(f, "Product_LRcorner_northing")
    ymax  <- hdf5r::h5attr(f, "Product_ULcorner_northing")
  }

  # get additional metadata ----
  sunzen  <- hdf5r::h5attr(f, "Sun_zenith_angle")
  sunaz   <- hdf5r::h5attr(f, "Sun_azimuth_angle")
  acqtime <- hdf5r::h5attr(f, "Product_StartTime")

  out_file_angles <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                            "_ANGLES.txt")
  utils::write.table(data.frame(date = acqtime,
                                sunzen   = sunzen,
                                sunaz = sunaz, stringsAsFactors = FALSE),
                     file = out_file_angles, row.names = FALSE)

  # get VNIR data cube and convert to raster ----

  out_file_vnir <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                          "_VNIR")
  out_file_vnir <- ifelse(out_format == "TIF",
                          paste0(out_file_vnir, ".tif"),
                          paste0(out_file_vnir, ".envi"))

  if (VNIR) {
    message("- Importing VNIR Cube -")
    if (file.exists(out_file_vnir) & !overwrite) {
      message("VNIR file already exists - use overwrite = TRUE or change output file name to reprocess")
      rast_vnir <- raster::stack(out_file_vnir)
    } else {

      if(proc_lev == 1) {
        vnir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source, "/Data Fields/VNIR_Cube")]][,,]
      } else {
        vnir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L2D_", source, "/Data Fields/VNIR_Cube")]][,,]
      }
      ind_vnir <- 1
      for (band_vnir in 1:66) {
        if (wl_vnir[band_vnir] != 0) {
          if(proc_lev == "1") {
            band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]), crs = "+proj=longlat +datum=WGS84")
            band <- raster::flip(band, 1)
            ex <- matrix(c(min(lon), max(lon),  min(lat), max(lat)), nrow = 2, ncol = 2, byrow = T)
            ex <- raster::extent(ex)
          } else {
            band <- raster::raster((vnir_cube[,order_vnir[band_vnir], ]),
                                   crs = paste0("+proj=utm +zone=", proj_code,
                                                " +datum=WGS84 +units=m +no_defs"))
            band <- t(band)
            ex <- matrix(c(xmin, xmax,  ymin, ymax), nrow = 2, ncol = 2, byrow = T)
            ex <- raster::extent(ex)
          }
          band <- raster::setExtent(band, ex, keepres=F)
          if (ind_vnir == 1) {
            rast_vnir <- band
          } else {
            rast_vnir <- raster::stack(rast_vnir, band)
          }
          ind_vnir <- ind_vnir + 1
        }
      }
      wl_vnir   <- wl_vnir[wl_vnir != 0]
      fwhm_vnir <- fwhm_vnir[fwhm_vnir != 0]
      names(rast_vnir) <- paste0("wl_", round(wl_vnir, digits = 4))
      rm(vnir_cube)
      rm(band)
      gc()

      if (VNIR) {
        message("- Writing VNIR raster -")
        rastwrite_lines(rast_vnir, out_file_vnir, out_format, proc_lev)
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

  out_file_swir <- paste0(tools::file_path_sans_ext(out_file),"_", source,
                          "_SWIR")
  out_file_swir <- ifelse(out_format == "TIF",
                          paste0(out_file_swir, ".tif"),
                          paste0(out_file_swir, ".envi"))

  if (SWIR) {
    if (file.exists(out_file_swir) & !overwrite) {
      message("SWIR file already exists - use overwrite = TRUE or change output file name to reprocess")
      rast_swir <- raster::stack(out_file_swir)
    } else {

      message("- Importing SWIR Cube - ")

      if (proc_lev == "1") {
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L1_", source, "/Data Fields/SWIR_Cube")]][,,]
      } else {
        swir_cube <- f[[paste0("HDFEOS/SWATHS/PRS_L2D_", source, "/Data Fields/SWIR_Cube")]][,,]
      }

      ind_band <- 1
      for (band_swir in 1:173) {
        if (wl_swir[band_swir] != 0) {
          if(proc_lev == "1") {
            band <- raster::raster((swir_cube[,order_swir[band_swir], ]), crs = "+proj=longlat +datum=WGS84")
            # band <- raster::t(raster::flip(band, 2))
            band <- raster::flip(band, 1)
            ex <- matrix(c(min(lon), max(lon),  min(lat), max(lat)), nrow = 2, ncol = 2, byrow = T)
            ex <- raster::extent(ex)
          } else {
            band <- raster::raster((swir_cube[,order_swir[band_swir], ]),
                                   crs = paste0("+proj=utm +zone=", proj_code,
                                                " +datum=WGS84 +units=m +no_defs"))
            band <- t(band)
            ex <- matrix(c(xmin, xmax,  ymin, ymax), nrow = 2, ncol = 2, byrow = T)
            ex <- raster::extent(ex)
          }
          band <- raster::setExtent(band, ex, keepres=F)

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

      if (SWIR) {
        message("- Writing SWIR raster -")
        rastwrite_lines(rast_swir, out_file_swir, out_format, proc_lev)
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
  out_file_full <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                          "_FULL")
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
        rast_tot <- raster::stack(rast_vnir[[which(wl_vnir < min(wl_swir))]], rast_swir)
        wl_tot   <- c(wl_vnir[which(wl_vnir < min(wl_swir))], wl_swir)
        fwhm_tot <- c(fwhm_vnir[which(wl_vnir < min(wl_swir))], fwhm_swir)
        names(rast_tot) <- paste0("wl_", round(wl_tot, digits = 4))
      } else {
        rast_tot <- raster::stack(rast_vnir, rast_swir[[which(wl_swir > max(wl_vnir))]])
        wl_tot   <- c(wl_vnir, wl_swir[which(wl_swir > max(wl_vnir))])
        fwhm_tot <- c(fwhm_vnir, fwhm_swir[which(wl_swir > max(wl_vnir))])
        names(rast_tot) <- paste0("wl_", round(wl_tot, digits = 4))
      }

      rm(rast_vnir)
      rm(rast_swir)
      gc()

      message("- Writing FULL raster -")
      rastwrite_lines(rast_tot, out_file_full, out_format, proc_lev)

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
  out_file_pan <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                         "_PAN")
  out_file_pan <- ifelse(out_format == "TIF",
                         paste0(out_file_pan, ".tif"),
                         paste0(out_file_pan, ".envi"))
  if (file.exists(out_file_pan) & !overwrite) {
    message("PAN file already exists - use overwrite = TRUE or change output file name to reprocess")
  } else {

    if (PAN) {

      message(" - Accessing PAN raster - ")

      if (proc_lev == "1") {
        pan_cube <- f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source), "/Data Fields/Cube")]][,]
        pan_lat <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                               "/Geolocation Fields/Latitude")]][,])
        pan_lon <- t(f[[paste0("/HDFEOS/SWATHS/PRS_L1_", gsub("H", "P", source),
                               "/Geolocation Fields/Longitude")]][,])
      } else {
        pan_cube <- f[[paste0("//HDFEOS/SWATHS/PRS_L2D_PCO/Data Fields/Cube")]][,]
      }
      if (proc_lev == "1") {
        rast_pan <- raster::raster(pan_cube, crs = "+proj=longlat +datum=WGS84")
        rast_pan <- raster::t(raster::flip(rast_pan, 2))
        ex   <- matrix(c(min(pan_lon), max(pan_lon),
                         min(pan_lat), max(pan_lat)),
                       nrow = 2, ncol = 2, byrow = T)
        ex   <- raster::extent(ex)
      } else {
        rast_pan <- raster::raster(pan_cube,
                                   crs = paste0("+proj=utm +zone=", proj_code,
                                                " +datum=WGS84 +units=m +no_defs"))
        ex <- matrix(c(xmin, xmax,  ymin, ymax), nrow = 2, ncol = 2, byrow = T)
        ex <- raster::extent(ex)
      }

      rm(pan_cube)
      gc()


      rast_pan <- raster::setExtent(rast_pan, ex, keepres=F)

      message("- Writing PAN raster -")

      rastwrite_lines(rast_pan, out_file_pan, out_format, proc_lev)
      rm(rast_pan)
      rm(pan_lon)
      rm(pan_lat)
      gc()
    }
  }

  if (proc_lev == 1) {
    # Save CLD if requested ----
    out_file_cld <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                           "_CLD")
    out_file_cld <- ifelse(out_format == "TIF",
                           paste0(out_file_cld, ".tif"),
                           paste0(out_file_cld, ".envi"))

    if (file.exists(out_file_cld) & !overwrite) {
      message("CLD file already exists - use overwrite = TRUE or change output file name to reprocess")
    } else {
      if (CLOUD) {

        message(" - Accessing CLOUD raster - ")
        cld_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/Cloud_Mask"]][,]
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

    # Save GLINT if requested ----
    out_file_glnt <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                            "_GLNT")
    out_file_glnt <- ifelse(out_format == "TIF",
                            paste0(out_file_glnt, ".tif"),
                            paste0(out_file_glnt, ".envi"))

    if (file.exists(out_file_glnt) & !overwrite) {
      message("GLINT file already exists - use overwrite = TRUE or change output file name to reprocess")
    } else {
      if (GLINT) {

        message(" - Accessing GLINT raster - ")
        glnt_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SunGlint_Mask"]][,]
        rast_glnt <- raster::raster(glnt_cube, crs = "+proj=longlat +datum=WGS84")
        rast_glnt <- raster::t(raster::flip(rast_glnt, 2))
        rm(glnt_cube)
        gc()

        ex   <- matrix(c(min(lon), max(lon),
                         min(lat), max(lat)),
                       nrow = 2, ncol = 2, byrow = T)
        ex   <- raster::extent(ex)
        rast_glnt <- raster::setExtent(rast_glnt, ex, keepres = F)

        message("- Writing GLINT raster -")
        rastwrite_lines(rast_glnt, out_file_glnt, out_format)

      }
    }

    # Save LC if requested ----
    out_file_lc <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                          "_LC")
    out_file_lc <- ifelse(out_format == "TIF",
                          paste0(out_file_lc, ".tif"),
                          paste0(out_file_lc, ".envi"))

    if (LC) {
      if (file.exists(out_file_lc) & !overwrite) {
        message("LC file already exists - use overwrite = TRUE or change output file name to reprocess")
      } else {

        message(" - Accessing LAND COVER raster - ")
        lc_cube <- f[["/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/LandCover_Mask"]][,]
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
}
