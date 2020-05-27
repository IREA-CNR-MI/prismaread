#' @title prisma_get_geoloc
#' @description helper function used to create accessory files needed for ATCOR
#'  atmospheric correction
#' @param wls passed by caller - array of PRISMA wavelengths
#' @param fwhms passed by caller - array of PRISMA fwhms
#' @param order_vnir passed by caller - ordering of array of PRISMA VNIR wavelengths
#' @param order_swir passed by caller - ordering of array of PRISMA SWIR wavelengths
#' @param f input data he5 from caller
#'
#' @inheritParams convert_prisma
#' @return `list` containing required info according to `proc_lev`
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#' @importFrom hdf5r existsGroup h5attr
prisma_make_atcor <- function(f,
                              out_file,
                              ATCOR_wls,
                              wls,
                              fwhms,
                              order_vnir,
                              order_swir,
                              join_priority,
                              source) {

    ATCOR_fold <- file.path(dirname(out_file), "ATCOR")
    out_file_wvl <- file.path(ATCOR_fold,
                              paste0(tools::file_path_sans_ext(basename(out_file)),
                                     "_atcor_wvl_nominal.wvl"))
    dir.create(ATCOR_fold, showWarnings = FALSE)

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

    if (!is.null(ATCOR_wls)) {

        if (source != "HCO") {

            if (hdf5r::existsGroup(f, "//KDP_AUX/Cw_Vnir_Matrix")) {

                if (!is.numeric(ATCOR_wls)) {
                    stop("ATCOR_wls should be either NULL or a vector containing the",
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
        } else {
            message("Creation of ATCOR files related to specific columns is only useful for HCR datasets - creation of additional ATCOR files skipped ")
        }
    }
}
