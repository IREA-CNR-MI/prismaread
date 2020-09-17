#' @title Compute spectral indexes from a PRISMA he5 file (or a file converted
#'  with pr_convert)
#' @description function used to compute spectral indexes, given the indexes
#'  formula
#' @details the function parses the index formula to identify the required
#'   bands. On the basis of identified bands, it retrieves the reflectance bands
#'   required, gets the data into R raster objects, performs the computation and
#'   stores results in a GeoTiff or ENVI raster file
#' @param in_file `character` path of the file to be used for computing indexes
#' @param out_file `character` output path - filenames are created adding a
#'  suffix to the basename of this file (e.g., S:/mypath/myoutfile_NDVI.tif)
#' @param out_format `character`` ["GTiff" | "ENVI"], Output format, Default:
#'  'GTiff'
#' @param indexes `character` array of names of indexes to be computed. You can
#'  see a list of available indexes using command `pr_listindexes()`, or see
#'  the corresponding table at:
#'  https://lbusett.github.io/prismaread/articles/Computing-Spectral-Indexes.html #nolint
#' @param cust_indexes `character` named list containing names and formulas of
#'  custom indexes to be computed. The indexes formulas must be computable R
#'  formulas, where bands are referred to by the prefix "b", followed by the
#'  wavelength (e.g., `cust_indexes = list(myindex1 = "R500 / R600",
#'                                  myindex2 = "(R800 - R680) / (R800 + R680)")`
#' @param overwrite `logical` if TRUE, existing files are overwritten,
#'  default: FALSE
#' @return NULL - new raster file saved in out_filename
#'
#' @author Lorenzo Busetto, phD (2017) \email{lbusett@@gmail.com}
#' @author Luigi Ranghetti, phD (2017) \email{ranghetti.l@@irea.cnr.it}
#' @note License: GPL 3.0
#' @export
#' @importFrom assertthat see_if
#' @importFrom tools file_path_sans_ext
#' @importFrom raster brick raster blockSize writeStart getValues writeValues
#'  writeStop
#' @importFrom stringr str_extract_all
#' @importFrom utils write.table
pr_compute_indexes <- function(in_file,
                               out_file,
                               out_format = "GTiff",
                               indexes,
                               cust_indexes = NULL,
                               overwrite    = FALSE) {
    if (!file.exists(in_file)) {
        stop("Input file not found. Aborting!")
    }
    assertthat::see_if(out_format %in% c("ENVI", "GTiff"))
    assertthat::see_if(dir.exists(dirname(out_file)),
                       msg = stop(
                           "Folder:", dirname(dirname(out_file)),
                           " does not exist. Please create it beforehand!"))

    out_file_txt <- paste0(tools::file_path_sans_ext(in_file), ".wvl")

    if (!file.exists(out_file_txt)) {
        stop("Unable to find the wavelengths txt file. Aborting!")
    }
    rast_wls <- read.table(out_file_txt, header = TRUE)$wl
    # create folder for index
    dir.create(dirname(out_file),
               showWarnings = FALSE, recursive = TRUE)

    in_rast <- raster::brick(in_file)

    index_list <- read.table(system.file("extdata/indexes_list.txt",
                                         package = "prismaread"), sep = "\t",
                             header = TRUE)
    av_indexes <- as.list(index_list$Formula)
    names(av_indexes) <- index_list$Name

    sel_indexes <- which(names(av_indexes) %in% indexes)

    bad_indexes <- which(!indexes %in% names(av_indexes))

    if (length(bad_indexes)) {
        warning("index(es) ", indexes[bad_indexes],
                " are not in the current list of ",
                "available indexes and will not be computed!")

    }

    tot_indexes <- c(av_indexes[sel_indexes], cust_indexes)
    indstrings <- list()
    for (ind in seq_along(tot_indexes)) {

        message("Computing ", names(tot_indexes)[[ind]], " Index")

        if (names(tot_indexes)[[ind]] %in% names(av_indexes)) {

            # create index filename
            out_indfile <- paste0(tools::file_path_sans_ext(out_file),
                                  "_", names(tot_indexes)[[ind]])
            out_indfile <- ifelse(out_format == "GTiff",
                                  paste0(out_indfile, ".tif"),
                                  paste0(out_indfile, ".envi"))
            # get index formula
            indform <- tot_indexes[[ind]]

            # compute index
            if (!file.exists(out_indfile) | overwrite) {
                req_wls <- as.numeric(substring(
                    stringr::str_extract_all(indform, "R[0-9,.]*")[[1]],2,100))
                which_bands <- unlist(lapply(req_wls,
                                             FUN = function(x) which.min(
                                                 abs(x - rast_wls))))

                diffs <- lapply(req_wls, FUN = function(x) min(
                    abs(x - rast_wls)))

                if (max(unlist(diffs)) > 9) {
                    warning("The bands in the input file do not allow to 2",
                            "compute index: ",
                            names(indexes)[[ind]])
                } else {
                    indstring    <- indform
                    in_rast_comp <- in_rast[[which_bands]]
                    rast_wls_comp <- rast_wls[which_bands]
                    which_bands  <- unlist(lapply(req_wls,
                                                  FUN = function(x) which.min(
                                                      abs(x - rast_wls_comp))))

                    for(nn in seq_along(req_wls)) {
                        indform   <- gsub(req_wls[nn],
                                          paste0("[,", which_bands[nn], "]"),
                                          indform)
                        indstring <- gsub(req_wls[nn],
                                          round(rast_wls_comp[nn], digits = 4),
                                          indstring)
                    }

                    indform <- gsub("R", "x", indform)
                    indstrings[[names(tot_indexes)[[ind]]]] <- indstring
                    out <- raster::raster(in_rast)
                    bs <-  raster::blockSize(out)
                    out <- raster::writeStart(out,
                                              filename = out_indfile,
                                              overwrite = TRUE,
                                              options = c("COMPRESS=LZW"),
                                              datatype = "FLT4S")


                    for (i in 1:bs$n) {
                        message("Writing Block: ", i, " of: ", bs$n)
                        x <- raster::getValues(in_rast_comp, row = bs$row[i],
                                               nrows = bs$nrows[i])
                        if (inherits(x, "numeric")) {
                            x <- x
                        } else {
                            x <- eval(parse(text = indform))
                        }
                        out <- raster::writeValues(out, x, bs$row[i])
                    }
                    out <- raster::writeStop(out)
                }
            } else {
                message("Output ", names(indexes)[[ind]],
                        " file already exists - use overwrite = TRUE or change
                    output file name to reprocess")
            }

        } else {
            warning("Index ", names(indexes)[[ind]],
                    " is not in the current list of available indexes and will",
                    " not be computed. ")
        }

    }
    indstrings_dt <- as.data.frame(indstrings)
    out_file_formulas <- paste0(tools::file_path_sans_ext(in_file),
                                ".formulas")
    utils::write.table(indstrings_dt,
                       file = out_file_formulas, row.names = FALSE)
}
