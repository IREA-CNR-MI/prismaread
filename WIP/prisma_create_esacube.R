prepare_prismacube <- function(in_l1file, in_l2cfile, in_wvl, out_folder, out_format = "GTIFF") {

    out_l1_folder <- file.path(out_folder, "L1")
    out_l2_folder <- file.path(out_folder, "L2")
    dir.create(out_l1_folder)
    dir.create(out_l2_folder)

    selbands_vnir <- in_wvl[which(in_wvl < 1000)]
    selbands_swir <- in_wvl[which(in_wvl >= 1000)]

    pr_convert(in_l1file, out_l1_folder, in_L2_file = in_l2cfile, selbands_vnir = selbands_vnir,
                   selbands_swir = selbands_swir, VNIR = FALSE, SWIR = FALSE, ANGLES = TRUE, LATLON = TRUE, LC = TRUE, CLOUD = TRUE,
                   PAN = FALSE, FULL = TRUE, out_format = "TIF", overwrite = F, fill_gaps = FALSE)

    pr_convert(in_l2cfile, out_l2_folder, in_L2_file = in_l2cfile, selbands_vnir = selbands_vnir,
                   selbands_swir = selbands_swir, VNIR = FALSE, SWIR = FALSE, ANGLES = FALSE, LATLON = FALSE, LC = FALSE,
                   PAN = FALSE, FULL = TRUE, out_format = "TIF", overwrite = F, fill_gaps = FALSE)

    stack <- raster::stack(
        file.path(out_l1_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_HCO_FULL.tif")),
        file.path(out_l2_folder, paste0(tools::file_path_sans_ext(basename(in_l2cfile)), "_HCO_FULL.tif")),
        file.path(out_l1_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_HCO_ANG.tif")),
        file.path(out_l1_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_HCO_LATLON.tif")),
        file.path(out_l1_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_HCO_LC.tif")),
        file.path(out_l1_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_HCO_CLD.tif")))

    centroid <- sf::st_coordinates(sf::st_centroid(sf::st_as_sfc(sf::st_bbox(stack))))
    # DTM <- raster::getData("SRTM", lon = centroid[1], lat=centroid[2])
    DTM <- elevatr::get_elev_raster(stack, z = 12)
    crpdtm <- raster::resample(DTM, stack)
    stack <- raster::stack(stack, crpdtm)

    wvls <- read.table(file.path(out_l1_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_HCO_FULL.wvl")),
                       header = T)

    stacknames <- c(paste0("Radiance_", round(wvls$wl, 5)),
                    paste0("Reflectance_", round(wvls$wl, 5)),
                    "ViewZen", "Relaz", "SunZen",
                    "Latitude", "Longitude",
                    "LandCover",
                    "Cloud_Mask",
                    "Elevation")
    names(stack) <- stacknames
    outstackname <- file.path(out_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_STACK.tif"))

    raster::writeRaster(stack, outstackname, options = c("COMPRESS=LZW"), overwrite = TRUE)


    outmetaname  <- file.path(out_folder, paste0(tools::file_path_sans_ext(basename(in_l1file)),  "_META.txt"))
    outmeta <- data.frame(band_n     = 1:raster::nlayers(stack),
                          band_name  = stacknames,
                          wavelength = as.numeric(c(wvls$wl, wvls$wl, rep("-", 8))),
                          fwhm       = as.numeric(c(wvls$fwhm, wvls$fwhm, rep("-", 8))),
                          units      = c(rep("Wm-2sr-1um-1", 8), rep("Reflectance", 8), rep("Deg", 5) ,"-", "-", "m"),
                          stringsAsFactors = FALSE)
    utils::write.table(outmeta,
                       file = outmetaname, row.names = FALSE)


}


# crpdtm <- raster::resample(DTM, stack)
#
# inlat <- raster::stack(file.path(out_l1_folder,paste0(tools::file_path_sans_ext(basename(in_l1file)), "_HCO_LATLON.tif")))[[1]]
# inlon <- raster::stack(file.path(out_l1_folder,paste0(tools::file_path_sans_ext(basename(in_l1file)), "_HCO_LATLON.tif")))[[2]]
#
# lat = inlat[1,]
# lon = inlon[,1]
#
# dimCross <- ncdim_def(name='cross_shore', units='m', longname='cross-shore coordinate', vals=data.distance )
#
# lon1 <- ncvar_def("latitude", "degrees_east", values(inlat))
# lat2 <- ncvar_def("longitude", "degrees_north", values(inlon))
#
# in_l1file  <- "/home/lb/projects/ESA_PRISMA_S5P_Covid/IMAGES/1_DATA/HDF_L1/PRS_L1_STD_OFFL_20200425103711_20200425103715_0001.he5"
# in_l2cfile <- "/home/lb/projects/ESA_PRISMA_S5P_Covid/IMAGES/1_DATA//HDF_L2C/PRS_L2C_STD_20200425103711_20200425103715_0001.he5"
# in_wvl     <- c(416,440,494,550, 670,772,865,2313)
# out_folder <- "/home/lb/projects/ESA_PRISMA_S5P_Covid/IMAGES/testcubes_0425"
# prepare_prismacube(in_l1file,in_l2cfile, in_wvl, out_folder )
# # # dir.create(out_folder)
# # #
# # # #
# in_l1file  <- "/home/lb/projects/ASI-PRISCAV/3_IMAGES/1_DATA/HDF5_L1/PRS_L1_STD_OFFL_20200714101655_20200714101659_0001.he5"
# in_l2cfile <- "/home/lb/projects/ASI-PRISCAV/3_IMAGES/1_DATA/HDF5_L2C/PRS_L2C_STD_20200714101655_20200714101659_0001.he5"
# in_wvl     <- c(416,440,494,550, 670,772,865,2313)
# out_folder <- "/home/lb/projects/ESA_PRISMA_S5P_Covid/IMAGES/testcubes_0714"
# # dir.create(out_folder)
# # #
# prepare_prismacube(in_l1file,in_l2cfile, in_wvl, out_folder )


