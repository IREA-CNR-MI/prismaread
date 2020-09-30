#' @title Extract data from a PRISMA he5 file and convert to GTiff or ENVI
#' @description Access a PRISMA HDF5 file and convert it to ENVI or GeoTiff
#'  format
#' @param in_file `character` full path of input HDF5 file
#' @param out_folder `character` full path to output folder
#' @param out_filebase `character` if "auto", output file names are based on
#'  filename of the input prisma file, with appropriate suffixes. If a different
#'  string is provided, it is used as prefix instead of the input file name,
#'  Default: "auto"
#' @param out_format `character`` ["GTiff" | "ENVI"], Output format, Default:
#' 'GTiff'
#' @param base_georef `logical` if TRUE, apply base georeferencing on L1, L2B/C
#'  data, using the "Georeference from input GLT" procedure explained here:
#'  https://www.harrisgeospatial.com/docs/backgroundgltbowtiecorrection.html,
#'  Default: TRUE
#' @param fill_gaps `logical` if TRUE, when georeferencing on L1, L2B/C data,
#'  substitute missing values with results of a 3x3 focal filter on the
#'  georeferenced data, Default: FALSE
#' @param source `character` ["HC0" | "HRC"], Considered Data Cub,
#'  Default: 'HCO'
#' @param VNIR `logical` if TRUE, create the VNIR image, Default: FALSE
#' @param SWIR `logical` if TRUE, create the SWIR image, Default: FALSE
#' @param FULL `logical` if TRUE, create a single multispectral image from
#'  VNIR and SWIR, Default: FALSE
#' @param join_priority `character` ["VNIR" | "SWIR"], spectrometer to consider
#'  when join_spectra = TRUE, Default: SWIR - ignored if join_spectra is FALSE.
#'  Default: "SWIR"
#' @param ATCOR  logical` if TRUE, create the text files required to run ATCOR,
#'  Default: FALSE;
#' @param ATCOR_wls  NULL, or `numeric` If NULL the only ATCOR wvl file
#'  created is the one containing Nominal wavelengths. If a numeric vector is
#'  provided, then one different wvl file is created for each selected COLUMN
#'  selected (e.g., if providing `ATCOR_wls = c(200, 800)`, then the wavelengths
#'  and FWHMs related to columns 200 and 800 are saved.)
#' @param PAN   `logical` if TRUE, also save PAN data, default: FALSE
#' @param CLOUD `logical` if TRUE, also save CLOUD MASK mask data,
#'  Default: FALSE (Ignored for L2 data).
#'  It is coded as:
#'    0 for not cloudy pixel
#'    1 for cloudy pixel
#'    10 = for not of all previous classification
#'    255 = error
#' @param GLINT `logical` if TRUE, also save GLINT mask data, default: FALSE
#' (Ignored for L2 data)
#'   It is coded as:
#'   0 for not sun glint
#'   1 for sun glint
#'   10 for not of all previous classification
#'   255 = error
#' @param LC `logical` if TRUE, also save the LAND COVER data, default: FALSE
#'  (Ignored for L2 data)
#'  It is coded as:
#'   0 for water pixel
#'   1 for snow pixel (and ice)
#'   2 for not-vegetated land pixel :bare soil)
#'   3 for crop and rangeland pixel
#'   4 for forest pixel
#'   5 for wetland pixel
#'   6 for not-vegetated land pixel :urban component
#'  10 for not of all previous classification
#'  255 for error
#' @param ANGLES if TRUE, also save the ACQUISITION ANGLES data. ANGLES data are
#'  saved as a 3-band raster. Band 1 contains "viewzen_ang", Band 2 contains
#'  "relaz_ang" and Band 3 contains "solzen_ang"), default: FALSE
#' @param LATLON if TRUE, also save the LATITUDE and LONGITUDE data. LATLON data
#'  are saved as a 2-band raster. Band 1 contains "Lat", Band 2 contains "Lon",
#'  Default: FALSE
#' @param ERR_MATRIX `logical` if TRUE, also save the SATURATION ERROR MATRIX
#' Data, Default: FALSE
#'  SATURATION ERROR MATRIX is coded as:
#'   0=pixel ok;
#'   1=DEFECTIVE PIXEL from KDP
#'   2= Pixel in saturation.
#'   3= Pixel with lower radiometric accuracy, due to coregistration effects.
#'   4= Pixel becomes NaN or Inf during processing.
#' @param overwrite `logical` if TRUE, existing files are overwritten,
#' Default: FALSE
#' @param apply_errmatrix  `logical` if TRUE, set pixels with values above 0 in
#'  the SATERR_MATRIX data cube (for L1) or ERRMATRIX cube (for L2) to NA,
#'  Default: FALSE
#' @param in_L2_file `character` full path of an L2B/C file to be used to
#'  extract georeferencing info and angles for a corresponding L1 file. If not
#'  NULL, and `in_file` is a L1 file, the LAT and LON fields used for bowtie
#'  georeferencing are taken from the L2 file instead than from the L1 file.
#'  The ANGLES data are also retrieved from the L2 file if requested.
#' @param selbands_vnir `numeric array` containing wavelengths (in nanometers)
#'  of bands that should be extracted from the VNIR data cube. If not NULL, only
#'  the bands with wavelengths closest to these values are extracted,
#'  Default: NULL
#' @param selbands_swir `numeric array` containing wavelengths (in nanometers)
#'  of bands that should be extracted from the SWIR data cube. If not NULL, only
#'  the bands with wavelengths closest to these values are extracted,
#'  Default: NULL
#' @param indexes `character` array of names of indexes to be computed. You can
#'  see a list of available indexes using command `pr_listindexes()`, or see
#'  the corresponding table at:
#'  https://lbusett.github.io/prismaread/articles/Computing-Spectral-Indexes.html, #nolint
#'  Default:NULL
#' @param cust_indexes `character` named list containing names and formulas of
#'  custom indexes to be computed. The indexes formulas must be computable R
#'  formulas, where bands are referred to by the prefix "b", followed by the
#'  wavelength (e.g.,`cust_indexes = list(myindex1 = "R500 / R600",
#'                                myindex2 = "(R800 - R680) / (R800 + R680)")`,
#'  Default:NULL
#' @param keep_index_cube `logical`, if TRUE, and spectral indexes were
#'  computed, keep in the output folder the hyperspectral cube used to compute
#'  the indexes (containing only wavelengths required to compute the desired
#'  indexes). The file takes a "_INDEXES" suffix. If FALSE, put the file in
#'  tempdir so that it is deleted automatically afterwards. Default: FALSE
#' @return The function is called for its side effects
#' @examples
#' \dontrun{
#'
#' # Here we download the example file from github using `piggyback` if not already
#' #  downloaded.
#' # WARNING ! This may need a long time (1GB file)
#
#'testfile_l1 <- file.path(
#' system.file("testdata/", package = "prismaread"),
#'             "PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.he5")
#'
#'if (!file.exists(testfile_l1)) {
#'   message("Downloading test data - This may need a long time!")
#'   piggyback::pb_download(
#'            "PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.zip",
#'            repo = "lbusett/prismaread",
#'            dest = file.path(
#'            system.file("", package = "prismaread"), "/testdata"))
#'   zipfile <- file.path(
#'                     system.file("testdata/", package = "prismaread"),
#'                     "PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.zip")
#'   unzip(zipfile, exdir = dirname(testfile_l1))
#'   unlink(zipfile)
#' }
#' out_folder = file.path(tempdir(), "prismaread_l1")
#'
#' # Save separately the VNIR and SWIR spetral cubes, in Tiff format. Also save
#' # LATLON, CLOUD and PAN
#'
#'  pr_convert(in_file    = testfile_l1,
#'             out_folder = out_folder,
#'             out_format = "GTiff",
#'             VNIR       = TRUE,
#'             SWIR       = TRUE,
#'             LATLON     = TRUE,
#'             PAN        = TRUE,
#'             CLOUD      = TRUE)
#'
#'  # Save the full cube, prioritizing the SWIR spectrometer and save in ENVI
#'  # format. Also save Cloud and LC data. Also, use a custom prefix for output
#'  # files
#'
#'  pr_convert(in_file       = testfile_l1,
#'             out_folder    = out_folder,
#'             out_filebase  = "prismaout_2020_05-14",
#'             out_format    = "ENVI",
#'             FULL          = TRUE,
#'             join_priority = "SWIR",
#'             LC            = TRUE,
#'             CLOUD         = TRUE,
#'             overwrite     = TRUE
#'             )
#'
#'  # Save a full image, prioritizing the VNIR spectrometer and save in TIF
#'  # format. Also create ATCOR files, with Nominal wavelengths
#'
#'  pr_convert(in_file       = in_file,
#'             out_folder    = out_folder,
#'             out_format    = "GTiff",
#'             FULL          = TRUE,
#'             join_priority = "VNIR",
#'             ATCOR         = TRUE
#'            )
#'
#' # Save the VNIR cube only, asociating a L2 file to be able to retrieve
#' # ANGLES data
#' testfile_l1 <- file.path(system.file("testdata/", package = "prismaread"),
#'                        "PRS_L2D_STD_20200524103704_20200524103708_0001.he5")
#' testfile_l2 <- file.path(system.file("testdata/", package = "prismaread"),
#'                         "PRS_L2C_STD_20200524103704_20200524103708_0001.he5")
#'
#'# Here we download the example file from github if not already downloaded.
#'# WARNING ! This may need a long time (1GB file)
#'
#'if (!file.exists(testfile_l2)){
#'  message("Downloading test data - This may need a long time!")
#'  piggyback::pb_download(
#'                  "PRS_L2C_STD_20200524103704_20200524103708_0001.zip",
#'                  repo = "lbusett/prismaread",
#'                  dest = file.path(
#'                        system.file("", package = "prismaread"), "/testdata"))
#'  zipfile <- file.path(system.file("testdata/", package = "prismaread"),
#'                       "PRS_L2C_STD_20200524103704_20200524103708_0001.zip")
#'  unzip(zipfile, exdir = dirname(testfile_l1))
#'  unlink(zipfile)
#'}
#'
#'out_folder = file.path(tempdir(), "prismaread_l2")
#'
#'pr_convert(in_file    = in_file,
#'           in_L2_file = in_L2_file,
#'           out_folder = out_folder,
#'           out_format = "ENVI",
#'           VNIR       = TRUE,
#'           ANGLES     = TRUE)
#' }
#' @rdname pr_convert
#' @export
#' @importFrom hdf5r H5File h5attr
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#' @importFrom raster stack
#'
pr_convert <- function(in_file,
                       out_folder,
                       out_filebase    = "auto",
                       out_format      = "ENVI",
                       base_georef     = TRUE,
                       fill_gaps       = FALSE,
                       VNIR            = FALSE,
                       SWIR            = FALSE,
                       FULL            = FALSE,
                       source          = "HCO",
                       join_priority   = "SWIR",
                       ATCOR           = FALSE,
                       ATCOR_wls       = NULL,
                       PAN             = FALSE,
                       CLOUD           = FALSE,
                       LC              = FALSE,
                       GLINT           = FALSE,
                       ANGLES          = FALSE,
                       LATLON          = FALSE,
                       ERR_MATRIX      = FALSE,
                       apply_errmatrix = FALSE,
                       overwrite       = FALSE,
                       in_L2_file      = NULL,
                       selbands_vnir   = NULL,
                       selbands_swir   = NULL,
                       indexes         = NULL,
                       cust_indexes    = NULL,
                       keep_index_cube = FALSE) {

  .pr_convert <- function(in_file, out_folder,
                          out_filebase, out_format, base_georef,
                          fill_gaps, VNIR,
                          SWIR, FULL,
                          source, join_priority,
                          ATCOR, ATCOR_wls,
                          PAN, CLOUD,
                          LC, GLINT,
                          ANGLES, LATLON,
                          ERR_MATRIX, apply_errmatrix,
                          overwrite,in_L2_file,
                          selbands_vnir, selbands_swir,
                          indexes,
                          cust_indexes,
                          keep_index_cube) {

    # Perform checks on inputs and return the hdf object

    f <- pr_check_inputs(in_file, out_folder,
                         out_filebase, out_format, base_georef,
                         fill_gaps, VNIR,
                         SWIR, FULL,
                         source, join_priority,
                         ATCOR, ATCOR_wls,
                         PAN, CLOUD,
                         LC, GLINT,
                         ANGLES, LATLON,
                         ERR_MATRIX, apply_errmatrix,
                         overwrite,in_L2_file,
                         selbands_vnir, selbands_swir,
                         indexes,
                         cust_indexes,
                         keep_index_cube)
    proc_lev <- hdf5r::h5attr(f, "Processing_Level")
    if (out_filebase == "auto") {
      out_filebase <- tools::file_path_sans_ext(basename(in_file))
    }
    out_file <- file.path(out_folder, out_filebase)

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


    # If indexes need to be computed, retrieve the list of VNIR and SWIR ----
    # wavelengths required for the computataion and automatically fill
    # the selbands_vnir and selbands_swir variables
    if (!is.null(indexes) | !is.null(cust_indexes)) {

      if (proc_lev %in% c("1", "2B")){

        warning("Spectral indexes are usually meant to be computed on ",
                "reflectance data. Proceed with caution!")

      }

      index_list <- read.table(system.file("extdata/indexes_list.txt",
                                           package = "prismaread"), sep = "\t",
                               header = TRUE)
      av_indexes <- as.list(index_list$Formula)
      names(av_indexes) <- index_list$Name

      sel_indexes <- which(names(av_indexes) %in% indexes)
      tot_indexes <- c(av_indexes[sel_indexes], cust_indexes)

      # when computing indexes, find out the required wavelengths
      # on vnir and swir
      getwls <- function(x) {
        substring(stringr::str_extract_all(x, "R[0-9,.]*")[[1]],2,100)
      }
      req_wls <- sort(as.numeric(
        unique(unlist((lapply(tot_indexes, FUN = function(x) getwls(x)))))))

      selbands_vnir <- req_wls[req_wls <= max(wl_vnir)]
      selbands_swir <- req_wls[req_wls >= min(wl_swir[wl_swir != 0])]

      if (any(selbands_swir %in% selbands_vnir)) {
        if (join_priority == "VNIR") {
          selbands_swir <- selbands_swir[selbands_swir >= max(wl_vnir)]
        } else {
          selbands_vnir <- selbands_vnir[selbands_vnir <=
                                           min(wl_swir[wl_swir != 0])]
        }
      }

      FULL <- TRUE
    }

    # write ATCOR files if needed ----
    if (ATCOR == TRUE && proc_lev == "1") {
      pr_make_atcor(f,
                    out_file,
                    ATCOR_wls,
                    wls,
                    fwhms,
                    order_vnir,
                    order_swir,
                    join_priority,
                    source)
    }


    # get additional metadata and create the "META" ancillary txt file----
    sunzen  <- hdf5r::h5attr(f, "Sun_zenith_angle")
    sunaz   <- hdf5r::h5attr(f, "Sun_azimuth_angle")
    acqtime <- hdf5r::h5attr(f, "Product_StartTime")

    out_file_angles <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                              ".ang")
    utils::write.table(data.frame(date = acqtime,
                                  sunzen   = sunzen,
                                  sunaz = sunaz, stringsAsFactors = FALSE),
                       file = out_file_angles, row.names = FALSE)


    # get VNIR data cube and convert to raster ----
    if (VNIR) {
      out_file_vnir <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                              "_VNIR")
      out_file_vnir <- ifelse(out_format == "GTiff",
                              paste0(out_file_vnir, ".tif"),
                              paste0(out_file_vnir, ".envi"))
    } else {
      out_file_vnir <- file.path(tempdir(), paste0(tools::file_path_sans_ext(
        basename(out_file)), "_", source,
        "_VNIR"))
      out_file_vnir <- paste0(out_file_vnir, ".tif")
    }

    if (VNIR | FULL) {
      message("- Importing VNIR Cube -")
      if (file.exists(out_file_vnir) & !overwrite) {
        message("VNIR file already exists - use overwrite = TRUE or change
                output file name to reprocess")
        rast_vnir <- raster::stack(out_file_vnir)
      } else {

        pr_create_vnir(f,
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
                       selbands_vnir = selbands_vnir,
                       in_L2_file = in_L2_file)
      }
    }

    # Build array of effectively processed bands/wl/fwhm  ----
    if (!is.null(selbands_vnir)){
      seqbands_vnir <- unlist(
        lapply(selbands_vnir, FUN = function(x) which.min(abs(wl_vnir - x))))

    } else {
      seqbands_vnir <- (1:66)[wl_vnir != 0]
    }
    wl_vnir   <- wl_vnir[seqbands_vnir]
    fwhm_vnir <- fwhm_vnir[seqbands_vnir]

    # get SWIR data cube and convert to raster ----

    if (SWIR) {
      out_file_swir <- paste0(tools::file_path_sans_ext(out_file),"_", source,
                              "_SWIR")
      out_file_swir <- ifelse(out_format == "GTiff",
                              paste0(out_file_swir, ".tif"),
                              paste0(out_file_swir, ".envi"))
    } else {
      out_file_swir <- file.path(tempdir(), paste0(tools::file_path_sans_ext(
        basename(out_file)), "_", source,
        "_SWIR"))
      out_file_swir <- paste0(out_file_swir, ".tif")
    }

    if ((SWIR | FULL) &&
        (is.null(selbands_swir) | length(selbands_swir) != 0)) {
      if (file.exists(out_file_swir) & !overwrite) {
        message("SWIR file already exists - use overwrite = TRUE or change
                output file name to reprocess")
        rast_swir <- raster::stack(out_file_swir)
      } else {

        message("- Importing SWIR Cube - ")
        pr_create_swir(f,
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
                       selbands_swir = selbands_swir,
                       in_L2_file = in_L2_file)
      }

      if (!is.null(selbands_swir)){
        seqbands_swir <- unlist(lapply(
          selbands_swir, FUN = function(x) which.min(abs(wl_swir - x))))
      } else {
        seqbands_swir <- (1:173)[wl_swir != 0]
      }
      wl_swir   <- wl_swir[seqbands_swir]
      fwhm_swir <- fwhm_swir[seqbands_swir]
    }

    # create FULL data cube and convert to raster ----
    if (FULL & is.null(indexes)) {
      out_file_full <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                              "_FULL")
      out_file_full <- ifelse(out_format == "GTiff",
                              paste0(out_file_full, ".tif"),
                              paste0(out_file_full, ".envi"))
    } else {
      if (keep_index_cube) {
        out_file_full <- paste0(tools::file_path_sans_ext(out_file), "_",
                                source,
                                "_INDEXES")
        out_file_full <- ifelse(out_format == "GTiff",
                                paste0(out_file_full, ".tif"),
                                paste0(out_file_full, ".envi"))
      } else {
        out_file_full <- file.path(tempdir(), paste0(tools::file_path_sans_ext(
          basename(out_file)), "_", source,
          "_INDEXES"))
        out_file_full <- paste0(out_file_full, ".tif")

      }
    }

    # Create and write the FULL hyperspectral cube if needed ----
    if (FULL) {
      if (file.exists(out_file_full) & !overwrite) {
        message("FULL file already exists - use overwrite = TRUE or change
                output file name to reprocess")
      } else {
        message("- Creating FULL raster -")

        # Save hyperspectral cube
        if (file.exists(out_file_vnir) && file.exists(out_file_swir)) {
          rast_vnir <- raster::stack(out_file_vnir)
          rast_swir <- raster::stack(out_file_swir)

          if (join_priority == "SWIR") {
            rast_tot <- raster::stack(
              rast_vnir[[which(wl_vnir < min(wl_swir))]], rast_swir)
            wl_tot   <- c(wl_vnir[which(wl_vnir < min(wl_swir))], wl_swir)
            fwhm_tot <- c(fwhm_vnir[which(wl_vnir < min(wl_swir))], fwhm_swir)
            names(rast_tot) <- c(
              paste0("b", seqbands_vnir[which(wl_vnir < min(wl_swir))], "_v"),
              paste0("b", seqbands_swir, "_s"))
          } else {
            rast_tot <- raster::stack(
              rast_vnir, rast_swir[[which(wl_swir > max(wl_vnir))]])
            wl_tot   <- c(wl_vnir, wl_swir[which(wl_swir > max(wl_vnir))])
            fwhm_tot <- c(fwhm_vnir, fwhm_swir[which(wl_swir > max(wl_vnir))])
            names(rast_tot) <- c(
              paste0("b", seqbands_vnir, "_v"),
              paste0("b", seqbands_swir[which(wl_swir > max(wl_vnir))], "_s"))
          }
          message("- Writing FULL raster -")

          pr_rastwrite_lines(rast_tot, out_file_full, out_format, proc_lev,
                             join = TRUE)
          rm(rast_vnir)
          rm(rast_swir)
          gc()
        } else {

          if (file.exists(out_file_vnir) && !file.exists(out_file_swir)) {
            rast_vnir <- raster::stack(out_file_vnir)
            message("SWIR file not created - FULL file will be equal to VNIR
                    one")
            file.copy(out_file_vnir, out_file_full)
            rast_tot <- rast_vnir
            wl_tot <- wl_vnir
            fwhm_tot <- fwhm_vnir
          } else {
            if (file.exists(out_file_swir) && !file.exists(out_file_vnir)) {
              rast_swir <- raster::stack(out_file_swir)
              message("VNIR file not created - FULL file will be equal to SWIR
                      one")
              file.copy(out_file_swir, out_file_full)
              rast_tot <- rast_swir
              wl_tot   <- wl_swir
              fwhm_tot <- fwhm_swir
            } else {
              warning("FULL file not created because neither SWIR nor VNIR
                      created")
            }
          }
        }

        # Write ENVI header if needed ----
        if (out_format == "ENVI") {
          cat("band names = {", paste(names(rast_tot),collapse=","), "}", "\n",
              file=raster::extension(out_file_full, "hdr"), append = TRUE)
          out_hdr <- paste0(tools::file_path_sans_ext(out_file_full), ".hdr")
          write(c("wavelength = {",
                  paste(round(wl_tot, digits = 4), collapse = ","), "}"),
                out_hdr, append = TRUE)
          write(c("fwhm = {",
                  paste(round(fwhm_tot, digits = 4), collapse = ","), "}"),
                out_hdr, append = TRUE)
          write("wavelength units = Nanometers", out_hdr, append = TRUE)
          write("sensor type = PRISMA", out_hdr, append = TRUE)
          write("data ignore value = -9.99000000e+002", out_hdr, append = TRUE)
        }

        out_file_txt <- paste0(tools::file_path_sans_ext(out_file_full), ".wvl")
        utils::write.table(data.frame(band   = seq_along(wl_tot),
                                      orband = substring(names(rast_tot) , 2),
                                      wl   = wl_tot,
                                      fwhm = fwhm_tot,
                                      stringsAsFactors = FALSE),
                           file = out_file_txt, row.names = FALSE)


        rm(rast_tot)
        gc()

      }

      # If only FULL selected, remove files related to VNIR and sWIR cubes if
      # existing
      if (!VNIR) {
        unlink(out_file_vnir)
      }

      if (!SWIR) {
        unlink(out_file_swir)
      }

      # now COMPUTE INDExEs if necessary ----
      if (!is.null(indexes) | !is.null(cust_indexes)) {

        pr_compute_indexes(in_file  = out_file_full,
                           out_file= out_file,
                           out_format = out_format,
                           indexes  = indexes,
                           cust_indexes = cust_indexes,
                           overwrite = overwrite)
      }

    }

    # Save PAN if requested ----
    out_file_pan <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                           "_PAN")
    out_file_pan <- ifelse(out_format == "GTiff",
                           paste0(out_file_pan, ".tif"),
                           paste0(out_file_pan, ".envi"))
    if (file.exists(out_file_pan) & !overwrite) {
      message("PAN file already exists - use overwrite = TRUE or change output
              file name to reprocess")
    } else {
      if (PAN) {
        pr_create_pan(f,
                      proc_lev,
                      source,
                      out_file_pan,
                      out_format,
                      base_georef,
                      fill_gaps,
                      in_L2_file = in_L2_file)

      }
    }

    # Save LATLON if requested ----
    out_file_latlon<- paste0(tools::file_path_sans_ext(out_file), "_", source,
                             "_LATLON")
    out_file_latlon <- ifelse(out_format == "GTiff",
                              paste0(out_file_latlon, ".tif"),
                              paste0(out_file_latlon, ".envi"))

    if (file.exists(out_file_latlon) & !overwrite) {
      message("LATLON file already exists - use overwrite = TRUE or change ",
              "output file name to reprocess")
    } else {
      if (LATLON) {
        pr_create_latlon(f,
                         proc_lev,
                         out_file_latlon,
                         out_format,
                         base_georef,
                         fill_gaps,
                         in_L2_file = in_L2_file)

      }
    }

    if (proc_lev %in% c("2B", "2C", "2D") | (proc_lev == "1" &
                                             !is.null(in_L2_file))) {
      # Save ANGLES if requested ----
      out_file_ang<- paste0(tools::file_path_sans_ext(out_file), "_", source,
                            "_ANG")
      out_file_ang <- ifelse(out_format == "GTiff",
                             paste0(out_file_ang, ".tif"),
                             paste0(out_file_ang, ".envi"))

      if (file.exists(out_file_ang) & !overwrite) {
        message("ANG file already exists - use overwrite = TRUE or change ",
                "output file name to reprocess")
      } else {
        if (ANGLES) {
          pr_create_angles(f,
                           proc_lev,
                           out_file_ang,
                           out_format,
                           base_georef,
                           fill_gaps,
                           in_L2_file = in_L2_file)

        }
      }
    }


    if (proc_lev == 1) {

      # Save CLD if requested ----
      out_file_cld <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                             "_CLD")
      out_file_cld <- ifelse(out_format == "GTiff",
                             paste0(out_file_cld, ".tif"),
                             paste0(out_file_cld, ".envi"))

      if (file.exists(out_file_cld) & !overwrite) {
        message("CLD file already exists - use overwrite = TRUE or change ",
                "output file name to reprocess")
      } else {
        if (CLOUD) {
          pr_create_additional(f,
                               type = "CLD",
                               out_file_cld,
                               out_format,
                               base_georef,
                               fill_gaps,
                               in_L2_file = in_L2_file)

        }
      }
      # Save GLINT if requested ----
      out_file_glnt <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                              "_GLINT")
      out_file_glnt <- ifelse(out_format == "GTiff",
                              paste0(out_file_glnt, ".tif"),
                              paste0(out_file_glnt, ".envi"))

      if (file.exists(out_file_glnt) & !overwrite) {
        message("GLINT file already exists - use overwrite = TRUE or change " ,
                "output file name to reprocess")
      } else {
        if (GLINT) {
          pr_create_additional(f,
                               type = "GLINT",
                               out_file_glnt,
                               out_format,
                               base_georef,
                               fill_gaps,
                               in_L2_file = in_L2_file)
        }
      }
      # Save LC if requested ----
      out_file_lc <- paste0(tools::file_path_sans_ext(out_file), "_", source,
                            "_LC")
      out_file_lc <- ifelse(out_format == "GTiff",
                            paste0(out_file_lc, ".tif"),
                            paste0(out_file_lc, ".envi"))

      if (LC) {
        if (file.exists(out_file_lc) & !overwrite) {
          message("LC file already exists - use overwrite = TRUE or change ",
                  "output file name to reprocess")
        } else {
          pr_create_additional(f,
                               type = "LC",
                               out_file_lc,
                               out_format,
                               base_georef,
                               fill_gaps,
                               in_L2_file = in_L2_file)
        }
      }
    }
  }

  # first run: ignore indexes ----

  .pr_convert(in_file,
              out_folder,
              out_filebase,
              out_format,
              base_georef,
              fill_gaps,
              VNIR,
              SWIR,
              FULL,
              source,
              join_priority,
              ATCOR,
              ATCOR_wls,
              PAN,
              CLOUD,
              LC,
              GLINT,
              ANGLES,
              LATLON,
              ERR_MATRIX,
              apply_errmatrix,
              overwrite,
              in_L2_file,
              selbands_vnir,
              selbands_swir,
              indexes       = NULL,
              cust_indexes  = NULL,
              keep_index_cube = FALSE)

  # second run: create indexes -----
  # in this way we can use the same function,
  # and when creating indexes create a temporary full raster, containing  only
  # bands required for the selected indexes

  if (!is.null(indexes)) {

    .pr_convert(in_file,
                out_folder,
                out_filebase,
                out_format,
                base_georef,
                fill_gaps,
                VNIR = FALSE,
                SWIR = FALSE,
                FULL = FALSE,
                source,
                join_priority,
                ATCOR = FALSE,
                ATCOR_wls,
                PAN = FALSE,
                CLOUD = FALSE,
                LC = FALSE,
                GLINT = FALSE,
                ANGLES = FALSE,
                LATLON = FALSE,
                ERR_MATRIX = FALSE,
                apply_errmatrix,
                overwrite,
                in_L2_file,
                selbands_vnir = NULL,
                selbands_swir = NULL,
                indexes       = indexes,
                cust_indexes  = cust_indexes,
                keep_index_cube = keep_index_cube)
  }


}
