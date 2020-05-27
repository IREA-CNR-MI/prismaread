#' @title convert_prisma
#' @description Access a PRISMA L1 HDF5 file and convert it to ENVI or GeoTiff
#'  format
#' @param in_file `character` full path of input HDF5 file
#' @param out_file `character` full path of output  file
#' @param out_format `character`` ["TIF" | "ENVI"], Output format, Default: 'tif'
#' @param base_georef `logical` if TRUE, apply base georeferencing on L1, L2B/C data,
#'  using the "Georeference from input GLT" procedure explained here:
#'  https://www.harrisgeospatial.com/docs/backgroundgltbowtiecorrection.html,
#'  Default: TRUE
#' @param fill_gaps `logical` if TRUE, when georeferencing on L1, L2B/C data,
#'  substitute missing values with results of a 3x3 focal filter on the georeferenced
#'  data, Default: TRUE
#' @param fix_geo `logical` if TRUE, apply a 90 metres correction to georeferencing
#' info of l2D data, Default: TRUE
#' @param source `character` ["HC0" | "HRC"], Considered Data Cube Default: 'HCO'
#' @param VNIR `logical` if TRUE, create the VNIR image, Default: TRUE
#' @param SWIR `logical` if TRUE, create the SWIR image, Default: TRUE
#' @param FULL `logical` if TRUE, create a single multispectral image from
#'  VNIR and SWIR, Default: FALSE
#' @param join_priority `character` ["VNIR" | "SWIR"], spectrometer to consider in
#'  the when join_spectra = TRUE, Default: SWIR - ignored if join_spectra is FALSE.
#'  Default: "VNIR"
#' @param ATCOR  logical` if TRUE, create the text files required to run ATCOR, Default: FALSE;
#' @param ATCOR_wls  NULL, or `numeric` If NULL the only ATCOR wvl file
#'  created is the one containing Nominal wavelengths. If a numeric vector is provided,
#'  then one different wvl file is created for each selected COLUMN selected (e.g., if
#'  providing `ATCOR_wls = c(200, 800)`, then the wavelengths and FWHMs related to
#'  columns 200 and 800 are saved.)
#' @param PAN   `logical` if TRUE, also save PAN data, default: TRUE (Ignored for L2 data)
#' @param CLOUD `logical` if TRUE, also save CLOUD MASK mask data, default: TRUE (Ignored for L2 data)
#' @param GLINT `logical` if TRUE, also save GLINT mask data, default: TRUE (Ignored for L2 data)
#' @param LC `logical` if TRUE, also save the LAND COVER data, default: TRUE (Ignored for L2 data)
#' @param ANGLES if TRUE, also save the ACQUISITION ANGLES data, default: TRUE (Ignored for L2 data)
#' @param LATLON if TRUE, also save the LATITUDA and LONGITUDE data, default: TRUE (Ignored for L2 data)
#' @param overwrite `logical` if TRUE, existing files are overwritten, default: FALSE
#' @param ERR_MATRIX Not yet implemented!
#' @param apply_errmatrix Not yet implemented!
#' @return The function is called for its side effects
#' @examples
#' \dontrun{
#'  in_file  <- "/home/lb/tmp/test/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
#'  out_file <- "/home/lb/tmp/test/test_1"
#'  out_format <- "ENVI"
#'
#'  # Save VNIR Cube image
#'  convert_prisma(in_file       = in_file,
#'                 out_file      = out_file,
#'                 out_format    = out_format,
#'                 VNIR          = TRUE,
#'                 SWiR          = TRUE,
#'                 FULL          = FALSE,
#'                 overwrite     = TRUE
#'                 )
#'
#'  # Save also the full cube, prioritizing the VNIR spectrometer and save in ENVI format
#'  # Also save Cloud and LC data
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
#' @importFrom utils write.table
#' @importFrom raster stack
#'
convert_prisma <- function(in_file,
                           out_folder,
                           out_filebase  = "auto",
                           out_format    = "ENVI",
                           base_georef   = TRUE,
                           fill_gaps     = TRUE,
                           fix_geo       = TRUE,
                           VNIR          = TRUE,
                           SWIR          = TRUE,
                           FULL          = FALSE,
                           source        = "HCO",
                           join_priority = "SWIR",
                           ATCOR         = FALSE,
                           ATCOR_wls     = NULL,
                           PAN           = TRUE,
                           CLOUD         = FALSE,
                           LC            = FALSE,
                           GLINT         = FALSE,
                           ANGLES        = FALSE,
                           LATLON        = FALSE,
                           ERR_MATRIX    = FALSE,
                           apply_errmatrix = FALSE,
                           overwrite     = FALSE) {

  # Open the file ----
  #
  if (!file.exists(in_file)) {
    stop("Selected input file does not exist. Verify your inputs. Aborting!")
  }

  f <- try(hdf5r::H5File$new(in_file, mode="r+"))

  if (inherits(f, "try-error")){
    stop("Unable to open the input file as a hdf5 file. Verify your inputs. Aborting!")
  }

  proc_lev <- hdf5r::h5attr(f, "Processing_Level")

  if (proc_lev != "1") {
    if (source %in% c("HRC", "HC0")) {
      message("Processing Level = 2 - Source modified to \"HCO\" by default")
      source = "HCO"
    }
  }

  if (!is.character(out_filebase)) {
    stop("out_filebase must be a string. Verify your inputs. Aborting!")
  }

  if (out_filebase == "auto") {
    out_filebase <- tools::file_path_sans_ext(basename(in_file))
  }

  out_file <- file.path(out_folder, out_filebase)

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
    prisma_make_atcor(f, out_file, ATCOR_wls, wls, fwhms, order_vnir, order_swir, join_priority, source)
  }

  # get geolocation info ----

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

  if (VNIR | FULL) {
    message("- Importing VNIR Cube -")
    if (file.exists(out_file_vnir) & !overwrite) {
      message("VNIR file already exists - use overwrite = TRUE or change output file name to reprocess")
      rast_vnir <- raster::stack(out_file_vnir)
    } else {
      prisma_create_vnir(f,
                         proc_lev,
                         source,
                         out_file_vnir,
                         out_format,
                         base_georef,
                         fill_gaps,
                         fix_geo,
                         wl_vnir,
                         order_vnir,
                         fwhm_vnir,
                         apply_errmatrix,
                         ERR_MATRIX)
    }
  }


  # be sure to remove zeroes also if VNIR already present to avoid  ----
  # errors on creation of FULL
  fwhm_vnir <- fwhm_vnir[wl_vnir != 0]
  wl_vnir <- wl_vnir[wl_vnir != 0]

  # get SWIR data cube and convert to raster ----

  out_file_swir <- paste0(tools::file_path_sans_ext(out_file),"_", source,
                          "_SWIR")
  out_file_swir <- ifelse(out_format == "TIF",
                          paste0(out_file_swir, ".tif"),
                          paste0(out_file_swir, ".envi"))

  if (SWIR | FULL) {
    if (file.exists(out_file_swir) & !overwrite) {
      message("SWIR file already exists - use overwrite = TRUE or change output file name to reprocess")
      rast_swir <- raster::stack(out_file_swir)
    } else {

      message("- Importing SWIR Cube - ")
      prisma_create_swir(f,
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
                         ERR_MATRIX)
    }
  }

  # be sure to remove zeroes also if swir already present to avoid  ----
  # errors on creation of FULL if recycling an existing cube from file
  fwhm_swir <- fwhm_swir[wl_swir != 0]
  wl_swir   <- wl_swir[wl_swir != 0]

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
      if (file.exists(out_file_vnir) && file.exists(out_file_swir)) {
        rast_vnir <- raster::stack(out_file_vnir)
        rast_swir <- raster::stack(out_file_swir)
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

        rastwrite_lines(rast_tot, out_file_full, out_format, proc_lev, join = TRUE)

        if (out_format == "ENVI") {
          out_hdr <- paste0(tools::file_path_sans_ext(out_file_full), ".hdr")
          write(c("wavelength = {",
                  paste(round(wl_tot, digits = 4), collapse = ","), "}"),
                out_hdr, append = TRUE)
          write(c("fwhm = {",
                  paste(round(fwhm_tot, digits = 4), collapse = ","), "}"),
                out_hdr, append = TRUE)
        }
        out_file_txt <- paste0(tools::file_path_sans_ext(out_file_full), "_wavelengths.txt")
        utils::write.table(data.frame(band = 1:length(wl_tot),
                                      wl   = wl_tot,
                                      fwhm = fwhm_tot, stringsAsFactors = FALSE),
                           file = out_file_txt, row.names = FALSE)

        rm(rast_tot)
        gc()
      } else {
        warning("Unable to create FULL data cube because VNIR or SWIR cubes were not created!")
      }
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
      prisma_create_pan(f,
                        proc_lev,
                        source,
                        out_file_pan,
                        out_format,
                        base_georef,
                        fill_gaps,
                        fix_geo)

    }
  }

  # Save LATLON if requested ----
  out_file_latlon<- paste0(tools::file_path_sans_ext(out_file), "_", source,
                           "_ANG")
  out_file_latlon <- ifelse(out_format == "TIF",
                            paste0(out_file_latlon, ".tif"),
                            paste0(out_file_latlon, ".envi"))

  if (file.exists(out_file_latlon) & !overwrite) {
    message("CLD file already exists - use overwrite = TRUE or change output file name to reprocess")
  } else {
    if (LATLON) {
      prisma_create_latlon(f,
                           proc_lev,
                           out_file_latlon,
                           out_format,
                           base_georef,
                           fill_gaps,
                           fix_geo)

    }
  }

  if (proc_lev %in% c("2B", "2C", "2D")) {
    # Save ANGLES if requested ----
    out_file_ang<- paste0(tools::file_path_sans_ext(out_file), "_", source,
                          "_ANG")
    out_file_ang <- ifelse(out_format == "TIF",
                           paste0(out_file_ang, ".tif"),
                           paste0(out_file_ang, ".envi"))

    if (file.exists(out_file_ang) & !overwrite) {
      message("CLD file already exists - use overwrite = TRUE or change output file name to reprocess")
    } else {
      if (ANGLES) {
        prisma_create_angles(f,
                             proc_lev,
                             out_file_ang,
                             out_format,
                             base_georef,
                             fill_gaps,
                             fix_geo)

      }
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
        prisma_create_additional(f,
                                 type = "CLD",
                                 out_file_cld,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 fix_geo)

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
        prisma_create_additional(f,
                                 type = "GLNT",
                                 out_file_glnt,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 fix_geo)
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
        prisma_create_additional(f,
                                 type = "LC",
                                 out_file_lc,
                                 out_format,
                                 base_georef,
                                 fill_gaps,
                                 fix_geo)
      }
    }
  }
}
