#' @title pr_check_inputs
#' @description helper function used to check consistency of pr_convert arguments
#' @inheritParams pr_convert
#' @return Upon success of all checks, the input hdf file, opened using
#'  `hdf5r::H5File`
#' @rdname pr_check_inputs
#' @importFrom assertthat see_if
#' @importFrom hdf5r H5File h5attr
pr_check_inputs <- function(in_file, out_folder,
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
    # Open the file ----

    # check consistency of logical arguments ----
    assertthat::see_if(is.logical(base_georef))
    assertthat::see_if(is.logical(fill_gaps))
    assertthat::see_if(is.logical(VNIR))
    assertthat::see_if(is.logical(SWIR))
    assertthat::see_if(is.logical(FULL))
    assertthat::see_if(is.logical(ATCOR))
    assertthat::see_if(is.logical(PAN))
    assertthat::see_if(is.logical(CLOUD))
    assertthat::see_if(is.logical(LC))
    assertthat::see_if(is.logical(GLINT))
    assertthat::see_if(is.logical(ANGLES))
    assertthat::see_if(is.logical(LATLON))
    assertthat::see_if(is.logical(ERR_MATRIX))
    assertthat::see_if(is.logical(apply_errmatrix))
    assertthat::see_if(is.logical(overwrite))
    assertthat::see_if(is.logical(keep_index_cube))

    assertthat::see_if(is.character(out_filebase))

    assertthat::see_if(out_format %in% c("ENVI", "GTiff"))

    assertthat::see_if(out_format %in% c("HCO", "HRC"))

    assertthat::see_if(join_priority %in% c("VNIR", "SWIR"))

    assertthat::see_if(is.null(ATCOR_wls) || is.numeric(ATCOR_wls))

    assertthat::see_if(is.numeric(selbands_vnir))
    assertthat::see_if(is.numeric(selbands_swir))

    assertthat::see_if(is.character(indexes))

    # check if cust_indexes is a list. If so, check the names
    assertthat::see_if(is.null(cust_indexes) || is.list(cust_indexes))
    if (!is.null(cust_indexes)){
        assertthat::see_if(is.null(pr_check_formulas(cust_indexes)) == TRUE)
    }

    f <- try(hdf5r::H5File$new(in_file, mode="r+"))

    if (inherits(f, "try-error")){
        stop("Unable to open the input file as a hdf5 file. Verify your inputs.
           Aborting!")
    }
    proc_lev <- hdf5r::h5attr(f, "Processing_Level")

    if (proc_lev != "1") {
        if (source %in% c("HRC")) {
            message("Processing Level = 2 - Source modified to \"HCO\" by ",
                    "default")
            source = "HCO"
        }
    }

    if (!is.character(out_filebase)) {
        stop("out_filebase must be a string. Verify your inputs. Aborting!")
    }

    out_file <- file.path(out_folder, out_filebase)

    if (!dir.exists(dirname(out_file))) {
        if (dir.exists(dirname(dirname(out_file)))) {
            dir.create(dirname(out_file))
        } else {
            stop("Folder:", dirname(dirname(out_file)),
                 " does not exist. Please create it beforehand!")
        }
    }
    if (!file.exists(in_file)) {
        stop("Selected input file does not exist. Verify your inputs. ",
             "Aborting!")
    }

    return(f)
}
