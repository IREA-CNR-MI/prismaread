#' @title pr_extract_spectra
#' @description function to extract values of a PRISMA image converted
#'  with prismaread on points/polygons saved on a vector file or on
#'  a `sf` object
#' @param in_file input PRISMA file obtained with `pr_convert`
#' @param in_vect either the full path to a vector file, or a `sf` object containing
#'  the points/polygons on which data should be extracted
#' @param id_field `character` (Optional), name of the column of the vector
#'  dataset to be used to "name" the outputs, and also "aggregate" them in case `dissolve` is TRUE.
#'  If NULL, a arbitrary `ID` field is created, and each point/polygon
#'  is considered separately (see Details),  Default: NULL
#' @param dissolve `logical` If TRUE and `id_field` was specified, in case multiple features of the input
#'  vector share a common id, they are dissolved before extracting the data, Default: FALSE
#' @param stats `logical` IF TRUE, compute standard statistics (mean, min, max, sd, variation coefficient)
#'  on the vector features, Default: TRUE
#' @param selstats `character` containing the statistics to be computed. Possible values are:
#'   "mean", "stdev","variance","coefficient_of_variation",
#'   "min","max"
#' @param stats_format `character` ["long" | "wide"] defines the format used for statistics output.
#'  If "long", the output has one column for the ID of the feature, and one column for each statistic.
#'  If "wide", the output has one column for each ID/statistic couple (e.g., mean_id_1, stdev_id_1, mean_id_2, etcetera)
#' @param quantiles `logical`, if TRUE, also compute quantiles on the vector features. Computed quantiles
#'  are set using the `percs` argument, Default: FALSE
#' @param percs `(sorted) numeric array [0,1]` defines which quantiles should be computed if
#'  `quantiles` is TRUE, Default: c(5,25,50,75,95)
#' @param allpix `logical` IF TRUE, also save the values for all pixels of the `in_vect`
#'  features in the `allpix` slot of the output list, Default: FALSE
#' @param out_file `character` full path of an output file where results should be stored, with
#'  extension. Valid extensions are ".csv", ".xls", ".xlsx" and ".RData". If NULL, output
#'  is not saved, Default: NULL
#' @return format of the output varies based on arguments `allpix` and `stats`
#'  1. If stats = TRUE and allpix = FALSE: a `tibble` containing extracted statistics,
#'     for each feature of the input and each wavelength. Format depends on `stat_sformat`;
#'  2. If stats = FALSE and allpix = TRUE:  a `tibble` containing extracted raster values,
#'     for each pixel of each feature of the input and each wavelength;
#'  3. If stats = TRUE and allpix = TRUE:  a `list` in which the `stats` slot contains
#'     statistics, and the `allpix` slot contains pixel values;
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  in_file <- "D:/prismareaetd/L2D/testL2D_HCO_VNIR.envi"
#'  in_vect <- "D:/prismaread/test/testpoints_l2d_polys.gpkg"
#'  # extract base statistics
#'  test <- pr_extract_spectra(in_file, in_vect, out_file = "D:/Temp/test1.xlsx")
#'  test
#'  # plot results using ggplot
#'  ggplot(test, aes(x = wvl, y = mean)) +
#'    geom_line(aes(color = ID, group = ID)) +
#'    facet_wrap(~ID) +
#'    theme_light()
#'
#'  # extract base statistics ands save results as excel file, in "wide" format
#'  test <- pr_extract_spectra(in_file, in_vect, out_file = "D:/Temp/test1.xlsx",
#'                                 stats_format = "wide")
#'  test
#'
#'  # extract custom statistics
#'  test <- pr_extract_spectra(in_file, in_vect,
#'                                 selstats = c("mean", "coeffvar", "stdev", "min", "max"))
#'  # plot results using ggplot
#'  ggplot(test, aes(x = wvl)) +
#'    geom_line(aes(y = mean, color = ID, group = ID)) +
#'    geom_line(aes(y = mean + stdev, group = ID), color = "grey75") +
#'    geom_line(aes(y = mean - stdev, group = ID), color = "grey75") +
#'    facet_wrap(~ID) +
#'    theme_light()
#'
#'  # extract custom statistics and quantiles
#'  test <- pr_extract_spectra(in_file, in_vect, quantiles = TRUE,
#'                                 selstats = c("mean", "stdev"))
#'  test
#'
#'  # extract also all pixels
#'  test <- pr_extract_spectra(in_file, in_vect, allpix = TRUE, quantiles = TRUE,
#'                                 selstats = c("mean", "stdev"))
#'  test$allpix
#'
#'  ggplot(test$allpix, aes(x = wvl)) +
#'    geom_line(aes(y = value, group = pixel, color = ID), lwd = 0.01)  +
#'    facet_wrap(~ID) +
#'    theme_light()
#'
#'  }
#' }
#' @rdname pr_extract_spectra
#' @export
#' @importFrom tools file_path_sans_ext file_ext
#' @importFrom raster brick res
#' @importFrom utils tail read.table write.csv
#' @importFrom data.table tstrsplit rbindlist
#' @importFrom sf st_read st_transform st_crs st_dimension st_buffer
#' @importFrom dplyr group_by summarise arrange filter select ungroup left_join
#' @importFrom rlang sym
#' @importFrom exactextractr exact_extract
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stringr str_pad
#' @importFrom tidyselect everything
#' @importFrom openxlsx write.xlsx
#' @importFrom stats quantile

pr_extract_spectra <- function(in_file,
                               in_vect,
                               id_field  = NULL,
                               dissolve  = FALSE,
                               stats     = TRUE,
                               selstats  = c("mean", "stdev"),
                               stats_format = "long",
                               quantiles = FALSE,
                               percs     =  c(0.05, 0.25,0.50,0.75,0.95),
                               allpix    = FALSE,
                               out_file  = NULL) {

  . <- pixel <- ID <- wvl <- value <- var <- NULL

  # Get the raster dataset ----
  if (!file.exists(in_file)) {
    stop("Input file does not exist. Aborting!")
  }

  if (!is.null(out_file)) {
    if(!dir.exists(dirname(out_file))) {
      stop("Folder specified for the output does not exist.
                 Please create it beforehand! Aborting!")
    } else {
      basefilename <- tools::file_path_sans_ext(out_file)
      ext <- tools::file_ext(out_file)
      if (!(ext %in% c("RData", "csv", "xls", "xlsx"))){
        stop("Extension for output file must be
                     \"RData\", \"csv\", \"xls\" or \"xlsx\". Aborting!")
      }
    }
  }

  if (!is.null(out_file)) {
    if(!dir.exists(dirname(out_file))) {
      stop("Folder specified for the output does not exist.
                 Please create it beforehand! Aborting!")
    } else {
      basefilename <- tools::file_path_sans_ext(out_file)
      ext <- tools::file_ext(out_file)
    }
  }

  in_rast <- try(raster::brick(in_file))
  in_type <- utils::tail(strsplit(basename(tools::file_path_sans_ext(in_file)),
                                  "_",
                                  fixed = TRUE)[[1]], 1)

  if (!in_type %in% c("VNIR", "SWIR", "FULL", "PAN", "LC", "CLD", "GLNT",
                      "ANGLES", "LATLON")) {
    stop("Input file does not seem to be a PRISMA file obtained from ",
         "PRISMAREAD. Aborting!")
  }

  if (in_type %in% c("VNIR", "SWIR", "FULL")) {
    in_file_wvl <- paste0(tools::file_path_sans_ext(in_file), ".wvl")
    wvl_ok <- FALSE
    if (tools::file_ext(in_file) == "envi"){
      # attempt to retrieve wavelengths from band names (works for ENVI files)
      tmpnames <- names(in_rast)
      tmpnames <- substring(tmpnames, 1, nchar(tmpnames) - 1)
      wvls <- as.numeric(data.table::tstrsplit(tmpnames, "..",
                                               fixed = TRUE)[[2]])
      if (is.numeric(wvls) && all(!is.na(wvls))) {
        wvl_ok <- TRUE
      } else {
        message("Unable to retrieve wavelengths from ENVI band names - ",
                "trying to use .wvl file")
      }
    }

    if (!wvl_ok) {
      if (file.exists(in_file_wvl)) {
        wvls <- utils::read.table(in_file_wvl, header = TRUE)
        wvls <- wvls$wl
        if (is.numeric(wvls) && all(!is.na(wvls))) {
          wvl_ok <- TRUE
        } else {
          message("Unable to retrieve wavelengths from .wvl file. Wavelengths set to
                            band index.")
          wvl <- 1:dim(in_rast)[3]
        }
      } else {
        message("Input Wavelengths file ", basename(in_file_wvl), " not found.",
                "Wavelengths set to band index. If you created the input file ",
                "with an older version of prismaread and want to use this ",
                "function, re-extract it to be able to retrieve wavelengths!")
      }
    }



  }

  if (!all(selstats %in% c("mean", "stdev", "variance", "min", "max",
                           "coeffvar"))) {
    stop("Invalid statistics requested. Supported statistics are:",
         "\"mean\", \"stdev\", \"variance\", \"min\", \"max\", \"coeffvar\". ",
         "Aborting!")
  }

  # Get the vector dataset ----
  if(is.character(in_vect)){
    if(file.exists(in_vect)) {
      in_sf <- try(sf::st_read(in_vect, quiet = TRUE))
      if(inherits(in_sf, "try-error")) {
        stop("in_vect does not appear to be a valid vector file. Aborting!")
      }
    } else {
      stop("in_vect file does not exist. Aborting!")
    }
  } else {
    if (!inherits(in_vect, "sf")) {
      stop("in_vect must be a `sf` object or a vector file in a GDAL-readable ",
           "format. Aborting!")
    } else {
      in_sf <- in_vect
    }
  }

  # dissolve features if needed ----
  if (dissolve && !is.null(id_field)){
    if(length(in_sf[[id_field]]) != length(unique(in_sf[[id_field]]))){
      in_sf <- in_sf %>%
        dplyr::group_by(!!rlang::sym(id_field))%>%
        dplyr::summarise()
    }
  }

  # extract the stats if needed ----
  in_sf <- sf::st_transform(in_sf, sf::st_crs(in_rast))

  # workaround to allow extraction over points
  if (all(sf::st_dimension(in_sf) == 0)) {
    in_sf <- sf::st_buffer(in_sf, raster::res(in_rast)[1] / 10000)
  }

  if (stats) {
    message("Extracting statistics")

    selstats_tmp <- selstats
    if ("coeffvar" %in% selstats) {
      selstats_tmp[which(selstats_tmp == "coeffvar")] <-
        "coefficient_of_variation"
    }
    out_vect <- t(exactextractr::exact_extract(in_rast, in_sf, selstats_tmp,
                                               progress = FALSE))

    if (!is.null(id_field)) {
      colnames(out_vect) <- in_sf[[id_field]]
    } else {
      colnames(out_vect) <- paste0("id_", seq_len(dim(in_sf)[1]))
    }



    varnames <- rep(selstats, each = dim(in_rast)[3])

    out_df <- data.frame(wvl = wvls, var = varnames, out_vect, row.names = NULL)

    out_df_w_stats <- out_df %>%
      tidyr::pivot_wider(., names_from = var, values_from =  3:dim(.)[2])

    out_df_l_stats <- out_df %>%
      tidyr::pivot_longer(., 3:dim(.)[2], names_to = "ID") %>%
      tidyr::pivot_wider(., names_from = c("var"), values_from = "value") %>%
      dplyr::arrange(ID, wvl)
  }
  # extract the pixels if needed----
  if (allpix | (stats & quantiles)) {
    message("Extracting pixel data")
    out_all_tmp <- exactextractr::exact_extract(in_rast, in_sf,
                                                progress = FALSE)
    names(out_all_tmp) <-   if (!is.null(id_field)) {
      in_sf[[id_field]]
    } else {
      paste0("id_", seq_len(dim(in_sf)[1]))
    }

    out_all   <- list()
    out_all_w <- list()
    for (ind in seq_along(out_all_tmp)) {

      tmp <- t(out_all_tmp[[ind]])
      out_df_all <- data.frame(wvl = c(wvls, "cov_frac"),
                               var = "value",
                               tmp, row.names = NULL)
      out_df_all_cov <- dplyr::filter(out_df_all, wvl == "cov_frac")
      out_df_w_all <- out_df_all %>%
        dplyr::filter(., wvl != "cov_frac") %>%
        tidyr::pivot_wider(., names_from = var, values_from =  3:dim(.)[2])
      colnames(out_df_w_all) <- c(
        "wvl",
        paste0("pix_",
               stringr::str_pad(seq_along(colnames(tmp)), max(nchar(colnames(tmp))),
                                "left", "0")))
      out_df_w_all$wvl <- as.numeric(as.character(out_df_w_all$wvl))
      # out_df_l_all <- out_df_all %>%
      #     tidyr::pivot_longer(., 3:dim(.)[2], names_to = "PIX")
      out_df_w_all$ID <- ifelse(is.null("id_field"), in_sf$id_field[[ind]],
                                paste0("id_", ind))
      out_df_w_all <- dplyr::select(out_df_w_all, ID, tidyselect::everything())
      #
      out_all_w[[ind]] <- out_df_w_all
      out_all[[ind]]   <- out_df_all
    }
    out_all_w <- data.table::rbindlist(out_all_w, fill = TRUE)
    out_all_l <- out_all_w %>%
      tidyr::pivot_longer(., 3:dim(.)[2], names_to = "pixel") %>%
      dplyr::arrange(ID, pixel, wvl)

    if (stats && quantiles) {
      for (qq in percs) {
        message("Computing quantiles - ", qq)
        perccol <- out_all_l %>%
          dplyr::group_by(ID, wvl) %>%
          dplyr::summarise(., perc = stats::quantile(value, probs  = qq,
                                              na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(., ID, wvl)
        names(perccol)[3] <- paste0("quant_", 100*qq)
        out_df_l_stats <- out_df_l_stats %>%
          dplyr::left_join(perccol, by = c("ID", "wvl"))
      }
    }
  } else {
    out_all_l = NULL
  }
  # Save if requested ----


  if(!allpix) {
    out_all_l <- NULL
  }

  if (stats && stats_format == "wide") {
    n_stats <- length(selstats)
    if (quantiles) {
      n_stats <- n_stats + length(percs)
    }
    n_ids <- length(unique(out_df_l_stats$ID))
    ordcol <- NULL
    for (st in 1:n_ids) {
      ordcol <- c(ordcol, 1 + st + (0:(n_stats-1)) * n_ids)
    }
    out_df_l_stats <- out_df_l_stats %>%
      tidyr::pivot_wider(., names_from = ID, values_from =  3:dim(.)[2]) %>%
      dplyr::select(c(1, ordcol))
  }
  # return ----

  if (stats & allpix) {
    out <- list("stats"  = out_df_l_stats,
                "allpix" = out_all_l)
  } else {
    if(!allpix) {
      out <- out_df_l_stats
    } else {
      out <- out_all_l
    }
  }

  # write to output file ----
  #
  if (!is.null(out_file)) {

    if (stats) {
      outfile_stats <- file.path(paste0(basefilename, "_stats.", ext))

      if (ext == "csv") {
        utils::write.csv(out_df_l_stats, file = outfile_stats,
                         row.names = FALSE)
      }
      if (ext %in% c("xls", "xlsx")) {
        openxlsx::write.xlsx(out_df_l_stats, file = outfile_stats,
                             rowNames = FALSE)
      }
      if (ext == "RData") {
        save(out_df_l_stats, file = outfile_stats)
      }
    }

    if (allpix) {
      outfile_all   <- file.path(paste0(basefilename, "_allpix.", ext))
      if (ext == "csv") {
        utils::write.csv(out_all_l, file = outfile_all, row.names = FALSE)
      }
      if (ext %in% c("xls", "xlsx")) {
        openxlsx::write.xlsx(out_all_l, file = outfile_all,
                             rowNames = FALSE)
      }
      if (ext == "RData") {
        save(out_all_l, file = outfile_all)
      }
    }
  }
  return(out)
}
