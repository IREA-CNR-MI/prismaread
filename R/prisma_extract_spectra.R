prisma_extract_spectra <- function(in_file,
                                   in_vect,
                                   id_field  = NULL,
                                   stats     = TRUE,
                                   quantiles = FALSE,
                                   allpix    = FALSE,
                                   out_file  = NULL) {

    # Get the raster dataset ----
    if (!file.exists(in_file)) {
        stop("Input file does not exist. Aborting!")
    }

    in_rast <- try(raster::brick(in_file))
    in_type <- tail(strsplit(basename(tools::file_path_sans_ext(in_file)), "_", fixed = TRUE)[[1]], 1)

    if (!in_type %in% c("VNIR", "SWIR", "FULL", "PAN", "LC", "CLD", "GLNT", "ANGLES", "LATLON")) {
        stop("Input file does not seem to be a PRISMA file obtained from PRISMAREAD. Aborting!")
    }
    if (!in_type %in% c("VNIR", "SWIR", "FULL")) {

        in_file_wvl <- paste0(tools::file_path_sans_ext(in_file), c("_wavelengths.txt", "_meta.txt"))
        if (any(file.exists(in_file_wvl))) {
            wvls <- read.table(in_file_wvl[which(file.exists(in_file_wvl))], header = TRUE)
        } else {
            stop("Input Wavelengths file ", basename(in_file_wvl), " not found. Aborting!")
        }
    }

    # Get the vector dataset ----
    if(is.character(in_vect)){
        if(file.exists(in_vect)) {
            in_sf <- try(sf::st_read(in_vect))
            if(inherits(in_sf, "try-error")) {
                stop("in_vect does not appear to be a valid vector file. Aborting!")
            }
        } else {
            stop("in_vect file does not exist. Aborting!")
        }
    } else {
        if (!inherits(in_vect, "sf")) {
            stop("in_vect must be a `sf` object or a vector file in a GDAL-readable format. Aborting!")
        }
    }

    # extract the stats if needed ----
    in_sf <- sf::st_transform(in_sf, sf::st_crs(in_rast))

    out_vect <- t(exactextractr::exact_extract(in_rast, in_sf, c("mean", "min", "max")))
    if (!is.null(id_field)) {
        colnames(out_vect) <- in_sf[[id_field]]
    } else {
        colnames(out_vect) <- paste0("id_", seq_len(dim(in_sf)[1]))
    }
    varnames <- rep(c("mean", "min", "max"), each = dim(in_rast)[3])
    wvls   <- wvls[,2]
    out_df <- data.frame(wvl = wvls, var = varnames, out_vect, row.names = NULL)

    out_df_w_stats <- out_df %>%
        tidyr::pivot_wider(., names_from = var, values_from =  3:dim(.)[2])

    out_df_l_stats <- out_df %>%
        tidyr::pivot_longer(., 3:dim(.)[2], names_to = "ID") %>%
        tidyr::pivot_wider(., names_from = c("var"), values_from = "value")

    # extract the pixels if needed----
    out_all_tmp <- exactextractr::exact_extract(in_rast, in_sf)
    names(out_all_tmp) <-   if (!is.null(id_field)) {
        in_sf[[id_field]]
    } else {
        paste0("id_", seq_len(dim(in_sf)[1]))
    }

    out_all   <- list()
    out_all_w <- list()
    for (ind in seq_along(out_all_tmp)) {
        tmp <- t(out_all_tmp[[ind]])
        colnames(tmp) <- paste0("pix_", seq_along(colnames(tmp)))
        out_df_all <- data.frame(wvl = c(wvls, "cov_frac"),  var = "value", tmp, row.names = NULL)
        out_df_all_cov <- dplyr::filter(out_df_all, wvl == "cov_frac")

        out_df_w_all <- out_df_all %>%
            dplyr::filter(., wvl != "cov_frac") %>%
            tidyr::pivot_wider(., names_from = var, values_from =  3:dim(.)[2])
        colnames(out_df_w_all) <- c("wvl", paste0("pix_", seq_along(colnames(tmp))))

        # out_df_l_all <- out_df_all %>%
        #     tidyr::pivot_longer(., 3:dim(.)[2], names_to = "PIX")
        out_df_w_all$ID <- ifelse(is.null("id_field"), in_sf$id_field[[ind]], paste0("id_", ind))
        out_df_w_all <- dplyr::select(out_df_w_all, ID, everything())
        #
        out_all_w[[ind]] <- out_df_w_all
        out_all[[ind]]   <- out_df_all
    }
    out_all_w <- data.table::rbindlist(out_all_w, fill = TRUE)
    out_all_l <- out_all_w %>%
        tidyr::pivot_longer(., 3:dim(.)[2], names_to = "pix")

    # Save if requested ----

}
