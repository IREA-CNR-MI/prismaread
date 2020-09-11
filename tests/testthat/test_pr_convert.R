# Tests on l2d -----

context("Access L2D data")
test_that(
    "Tests on L2D", {
        skip_on_cran()
        skip_on_travis()

        testfile <- file.path(system.file("testdata/", package = "prismaread"),
                              "PRS_L2D_STD_20200524103704_20200524103708_0001.he5")

        # Download and unzip using piggyback if necessary
        if (!file.exists(testfile)){
            message("Downloading test data - This may need a long time!")
            piggyback::pb_download("PRS_L2D_STD_20200524103704_20200524103708_0001.zip",
                                   repo = "lbusett/prismaread",
                                   dest = file.path(system.file("", package = "prismaread"), "/testdata"))
            piggyback::pb_track(glob = "inst/testdata/*.zip, inst/testdata/*.he5")
            zipfile <- file.path(system.file("testdata/", package = "prismaread"),
                                 "PRS_L2D_STD_20200524103704_20200524103708_0001.zip")
            unzip(zipfile, exdir = dirname(testfile))
            unlink(zipfile)
        }
        out_folder_L2D <- file.path(tempdir(), "prismaread/L2D")
        dir.create(out_folder_L2D, recursive = TRUE)

        # launch pr_convert to creat VNIR (FULL), SWIR (FULL), FULL,
        # ANGLES, LATLON, PAN - save to tiff
        pr_convert(in_file = testfile, out_folder = out_folder_L2D,
                   out_filebase = "testL2D_1", out_format = "GTiff",
                   VNIR = TRUE,
                   SWIR = TRUE, FULL = FALSE, ANGLES = FALSE, PAN = FALSE,
                   LATLON = FALSE, ERR_MATRIX = FALSE, apply_errmatrix = TRUE,
                   overwrite = TRUE)
        vnir  <- raster::brick(file.path(out_folder_L2D, "testL2D_1_HCO_VNIR.tif"))
        testthat::expect_equal(dim(vnir), c(1236,1290,63))
        swir  <- raster::brick(file.path(out_folder_L2D, "testL2D_1_HCO_SWIR.tif"))
        testthat::expect_equal(dim(swir), c(1236,1290,171))


        # launch pr_convert to creat VNIR (6 bands), SWIR (3 bands), FULL,
        # ANGLES, LATLON, PAN - save to tiff
        pr_convert(in_file = testfile, out_folder = out_folder_L2D,
                   out_filebase = "testL2D_1", out_format = "GTiff",
                   VNIR = TRUE, selbands_vnir = c(450, 550, 650, 750, 850, 950),
                   selbands_swir = c(1000, 1500 ,2300),
                   SWIR = TRUE, FULL = TRUE, ANGLES = TRUE, PAN = TRUE,
                   LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                   overwrite = TRUE)

        # All files created
        flist <- list.files(out_folder_L2D, pattern = "testL2D_1")
        testthat::expect_equal(length(flist), 12)

        vnir  <- raster::brick(file.path(out_folder_L2D, "testL2D_1_HCO_VNIR.tif"))
        means_vnir <- as.numeric(raster::cellStats(vnir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_vnir, c(0.022935, 0.05126269, 0.05515500,
                                             0.12456146, 0.14181177, 0.14477050))

        swir  <- raster::brick(file.path(out_folder_L2D, "testL2D_1_HCO_SWIR.tif"))
        means_swir <- as.numeric(raster::cellStats(swir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_swir, c(0.13954069, 0.10219114, 0.08076247))

        full <- raster::brick(file.path(out_folder_L2D, "testL2D_1_HCO_FULL.tif"))
        means_full <- as.numeric(raster::cellStats(full, "mean", na.rm = TRUE))
        testthat::expect_equal(means_full, c(0.022935, 0.05126269, 0.05515500,
                                             0.12456146, 0.14181177, 0.14477050,
                                             0.13954069, 0.10219114, 0.08076247))

        testthat::expect_true(raster::projection(vnir) %in% c(
            "+proj=longlat +datum=WGS84 +no_defs",
            "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"))

        # wavelengths in wvl files are correct
        wvls_vnir <- read.table(file.path(out_folder_L2D, "testL2D_1_HCO_VNIR.wvl"),
                                header = TRUE)
        wvls_swir <- read.table(file.path(out_folder_L2D, "testL2D_1_HCO_SWIR.wvl"),
                                header = TRUE)
        wvls_full <- read.table(file.path(out_folder_L2D, "testL2D_1_HCO_FULL.wvl"),
                                header = TRUE)
        testthat::expect_equal(wvls_full$wl, c(wvls_vnir$wl, wvls_swir$wl))

        # extent is correct
        ext <- raster::extent(vnir)
        testthat::expect_equal(ext[c(1,3)], c(443902.5 - 15, 5006667.5 -15))

        pan <- raster::brick(file.path(out_folder_L2D, "testL2D_1_HCO_PAN.tif"))
        rastres <- raster::res(pan)
        testthat::expect_equal(rastres, c(5, 5))
        # launch pr_convert to creat VNIR (3 bands), SWIR (3 bands), FULL,
        # - save to envi
        pr_convert(in_file = testfile, out_folder = out_folder_L2D,
                   out_filebase = "testL2D_2", out_format = "ENVI",
                   VNIR = TRUE, selbands_vnir = c(450, 650, 850),
                   selbands_swir = c(1500),
                   SWIR = TRUE, FULL = TRUE, ANGLES = FALSE,
                   LATLON = FALSE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                   overwrite = TRUE)

        flist <- list.files(out_folder_L2D, pattern = "testL2D_2")
        testthat::expect_equal(length(flist), 19)
        vnir  <- raster::brick(file.path(out_folder_L2D, "testL2D_2_HCO_VNIR.envi"))
        means_vnir <- as.numeric(raster::cellStats(vnir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_vnir, c(0.022935, 0.05515500, 0.14181177))

        swir  <- raster::brick(file.path(out_folder_L2D, "testL2D_2_HCO_SWIR.envi"))
        means_swir <- as.numeric(raster::cellStats(swir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_swir, c(0.10219114))

        # see that apply_errmatrix does something
        pr_convert(in_file = testfile, out_folder = out_folder_L2D,
                   out_filebase = "testL2D_3", out_format = "GTiff",
                   VNIR = TRUE, selbands_vnir = c(450, 650),
                   SWIR = FALSE, FULL = FALSE, ANGLES = FALSE,
                   LATLON = FALSE, ERR_MATRIX = FALSE, apply_errmatrix = TRUE,
                   overwrite = TRUE)
        vnir <- raster::brick(file.path(out_folder_L2D, "testL2D_3_HCO_VNIR.tif"))
        means_vnir <- as.numeric(raster::cellStats(vnir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_vnir, c(0.02349694, 0.05520516))

        # wavelengths in wvl files are correct
        wvls_vnir <- read.table(file.path(out_folder_L2D, "testL2D_1_HCO_VNIR.wvl"),
                                header = TRUE)
        wvls_swir <- read.table(file.path(out_folder_L2D, "testL2D_1_HCO_SWIR.wvl"),
                                header = TRUE)
        wvls_full <- read.table(file.path(out_folder_L2D, "testL2D_1_HCO_FULL.wvl"),
                                header = TRUE)
        testthat::expect_equal(wvls_full$wl, c(wvls_vnir$wl, wvls_swir$wl))
    })

# Tests on l2c -----

context("Access L2C data")
test_that(
    "Tests on L2C", {
        skip_on_cran()
        skip_on_travis()

        testfile <- file.path(system.file("/testdata", package = "prismaread"),
                              "/PRS_L2C_STD_20200524103704_20200524103708_0001.he5")

        # Download and unzip using piggyback if necessary
        if (!file.exists(testfile)){
            message("Downloading test data - This may need a long time!")
            piggyback::pb_download("PRS_L2C_STD_20200524103704_20200524103708_0001.zip",
                                   repo = "lbusett/prismaread",
                                   dest = file.path(system.file("", package = "prismaread"), "/testdata"))
            piggyback::pb_track(glob = "inst/testdata/*.zip, inst/testdata/*.he5")
            zipfile <- file.path(system.file("/testdata", package = "prismaread"),
                                 "/PRS_L2C_STD_20200524103704_20200524103708_0001.zip")
            unzip(zipfile, exdir = dirname(testfile))
            unlink(zipfile)
        }
        out_folder_L2C <- file.path(tempdir(), "prismaread/L2C")
        dir.create(out_folder_L2C, recursive = TRUE)

        # launch pr_convert to creat VNIR (6 bands), SWIR (3 bands), FULL,
        # ANGLES, LATLON, PAN - save to tiff - with georef and gaps filling
        pr_convert(in_file = testfile, out_folder = out_folder_L2C,
                   out_filebase = "testL2C_1", out_format = "GTiff",
                   VNIR = TRUE, selbands_vnir = c(450, 550, 650, 750, 850, 950),
                   selbands_swir = c(1000, 1500 ,2300),
                   SWIR = TRUE, FULL = TRUE, ANGLES = TRUE, PAN = TRUE,
                   LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                   overwrite = TRUE)

        # All files created
        flist <- list.files(out_folder_L2C, pattern = "testL2C_1")
        testthat::expect_equal(length(flist), 12)

        vnir  <- raster::brick(file.path(out_folder_L2C, "testL2C_1_HCO_VNIR.tif"))
        means_vnir <- as.numeric(raster::cellStats(vnir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_vnir, c(0.03527613  , 0.07889447 , 0.08488756  ,
                                             0.19185325  , 0.21844948  , 0.22302736),
                               tolerance = 0.001)

        swir  <- raster::brick(file.path(out_folder_L2C, "testL2C_1_HCO_SWIR.tif"))
        means_swir <- as.numeric(raster::cellStats(swir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_swir, c(0.2149935  , 0.1573347  , 0.1243109),
                               tolerance = 0.001)

        full <- raster::brick(file.path(out_folder_L2C, "testL2C_1_HCO_FULL.tif"))
        means_full <- as.numeric(raster::cellStats(full, "mean", na.rm = TRUE))
        testthat::expect_equal(means_full, c(0.03527613, 0.07889447, 0.08488756,
                                             0.19185325, 0.21844948, 0.22302736,
                                             0.2149935,  0.1573347, 0.1243109),
                               tolerance = 0.001)

        angs  <- raster::brick(file.path(out_folder_L2C, "testL2C_1_HCO_ANG.tif"))
        means_angs <- as.numeric(raster::cellStats(angs, "mean", na.rm = TRUE))
        testthat::expect_equal(means_angs, c(11.73522, 129.13642, 26.24126),
                               tolerance = 0.001)

        llon  <- raster::brick(file.path(out_folder_L2C, "testL2C_1_HCO_LATLON.tif"))
        means_latlon <- as.numeric(raster::cellStats(llon, "mean", na.rm = TRUE))
        testthat::expect_equal(means_latlon, c(45.379759, 8.529825),
                               tolerance = 0.001)

        # wavelengths in wvl files are correct
        wvls_vnir <- read.table(file.path(out_folder_L2C, "testL2C_1_HCO_VNIR.wvl"),
                                header = TRUE)
        wvls_swir <- read.table(file.path(out_folder_L2C, "testL2C_1_HCO_SWIR.wvl"),
                                header = TRUE)
        wvls_full <- read.table(file.path(out_folder_L2C, "testL2C_1_HCO_FULL.wvl"),
                                header = TRUE)
        testthat::expect_equal(wvls_full$wl, c(wvls_vnir$wl, wvls_swir$wl))

        # extent is correct: check left and top and not corner due to rotation
        # introduced by GLT
        ext <- raster::extent(vnir)
        testthat::expect_equal(ext[c(1,4)], c(8.2876167, 45.544937), tolerance = 0.00001)

        # launch pr_convert to creat VNIR (3 bands), SWIR (3 bands), FULL,
        # - save to GTiff - no fill gaps
        pr_convert(in_file = testfile, out_folder = out_folder_L2C,
                   out_filebase = "testL2C_2", out_format = "GTiff",
                   fill_gaps = FALSE,
                   VNIR = TRUE, selbands_vnir = c(450, 650, 850),
                   selbands_swir = c(1500),
                   SWIR = TRUE, FULL = TRUE, ANGLES = FALSE,
                   LATLON = FALSE, ERR_MATRIX = FALSE, apply_errmatrix = FALSE,
                   overwrite = TRUE)

        flist <- list.files(out_folder_L2C, pattern = "testL2C_2")
        testthat::expect_equal(length(flist), 7)

        vnir  <- raster::brick(file.path(out_folder_L2C, "testL2C_2_HCO_VNIR.tif"))
        means_vnir <- as.numeric(raster::cellStats(vnir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_vnir, c(0.03529589   , 0.08492260  , 0.21841796),
                               tolerance = 0.001)

        swir  <- raster::brick(file.path(out_folder_L2C, "testL2C_2_HCO_SWIR.tif"))
        means_swir <- as.numeric(raster::cellStats(swir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_swir, c(0.1573508), tolerance = 0.001)

        full <- raster::brick(file.path(out_folder_L2C, "testL2C_2_HCO_FULL.tif"))
        means_full <- as.numeric(raster::cellStats(full, "mean", na.rm = TRUE))
        testthat::expect_equal(means_full, c(0.03529589, 0.08492260, 0.21841796,
                                             0.1573508), tolerance = 0.001)

        testthat::expect_equal(raster::projection(vnir),
                               "+proj=longlat +datum=WGS84 +no_defs")

        # launch pr_convert to creat VNIR (3 bands), SWIR (3 bands), FULL,
        # - save to GTiff - no georef
        testthat::expect_warning(pr_convert(in_file = testfile, out_folder = out_folder_L2C,
                                            out_filebase = "testL2C_3", out_format = "GTiff",
                                            fill_gaps = FALSE, base_georef = FALSE,
                                            VNIR = TRUE, selbands_vnir = c(450, 650, 850),
                                            selbands_swir = c(1500),
                                            SWIR = TRUE, FULL = TRUE, ANGLES = TRUE,
                                            LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                                            overwrite = TRUE))

        flist <- list.files(out_folder_L2C, pattern = "testL2C_3")
        testthat::expect_equal(length(flist), 11)
        vnir  <- raster::brick(file.path(out_folder_L2C, "testL2C_3_HCO_VNIR.tif"))
        testthat::expect_error(projection(vnir))
        testthat::expect_equal(dim(vnir), c(1000,1000,3))

    })


context("Access L1 data")
# Tests on l1 -----
test_that(
    "Tests on L1", {
        skip_on_cran()
        skip_on_travis()

        testfile <- file.path(system.file("testdata/", package = "prismaread"),
                              "PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.he5")

        # Download and unzip using piggyback if necessary
        if (!file.exists(testfile)){
            message("Downloading test data - This may need a long time!")
            piggyback::pb_download("PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.zip",
                                   repo = "lbusett/prismaread",
                                   dest = file.path(system.file("", package = "prismaread"), "/testdata"))
            piggyback::pb_track(glob = "inst/testdata/*.zip, inst/testdata/*.he5")
            zipfile <- file.path(
                system.file("testdata/", package = "prismaread"),
                "PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.zip")
            unzip(zipfile, exdir = dirname(testfile))
            unlink(zipfile)
        }
        out_folder_L1 <- file.path(tempdir(), "prismaread/L1")
        dir.create(out_folder_L1, recursive = TRUE)

        # launch pr_convert to creat VNIR (3 bands), SWIR (1 bands), FULL,
        # ANGLES, LATLON, PAN - save to tiff
        pr_convert(in_file = testfile, out_folder = out_folder_L1,
                   out_filebase = "testL1_1", out_format = "GTiff",
                   VNIR = TRUE, selbands_vnir = c(450, 650, 850),
                   selbands_swir = c(1500), ATCOR = TRUE,
                   SWIR = TRUE, FULL = TRUE, ANGLES = TRUE, PAN = TRUE,
                   CLOUD = TRUE, GLINT = TRUE, LC = TRUE,
                   LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                   overwrite = TRUE)

        # All files created
        flist <- list.files(out_folder_L1, pattern = "testL1_1")
        testthat::expect_equal(length(flist), 14)

        vnir  <- raster::brick(file.path(out_folder_L1, "testL1_1_HCO_VNIR.tif"))
        means_vnir <- as.numeric(raster::cellStats(vnir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_vnir, c(71.72868  , 45.73275  , 58.20885),
                               tolerance = 0.0001)

        swir  <- raster::brick(file.path(out_folder_L1, "testL1_1_HCO_SWIR.tif"))
        means_swir <- as.numeric(raster::cellStats(swir, "mean", na.rm = TRUE))
        testthat::expect_equal(means_swir, c(8.921804), tolerance = 0.0001)

        full <- raster::brick(file.path(out_folder_L1, "testL1_1_HCO_FULL.tif"))
        means_full <- as.numeric(raster::cellStats(full, "mean", na.rm = TRUE))
        testthat::expect_equal(means_full, c(71.72868, 45.73275, 58.20885,
                                             8.921804), tolerance = 0.0001)

        # Use L2file to georeference ----
        in_L2_file <- file.path(
            system.file("testdata/", package = "prismaread"),
            "PRS_L2C_STD_20200524103704_20200524103708_0001.he5")
        pr_convert(in_file = testfile, out_folder = out_folder_L1,
                   out_filebase = "testL1_2", out_format = "GTiff",
                   in_L2_file = in_L2_file,
                   VNIR = TRUE, selbands_vnir = c(450, 650, 850),
                   selbands_swir = c(1500),
                   SWIR = TRUE, FULL = TRUE, ANGLES = TRUE, PAN = TRUE,
                   CLOUD = TRUE, GLINT = TRUE, LC = TRUE,
                   LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                   overwrite = TRUE)

        angs  <- raster::brick(file.path(out_folder_L1, "testL1_2_HCO_ANG.tif"))
        means_angs <- as.numeric(raster::cellStats(angs, "mean", na.rm = TRUE))
        # angles retrieved by overwriting with L2 are equal to those of L2
        testthat::expect_equal(means_angs, c(11.73522, 129.13642, 26.24126),
                               tolerance = 0.001)

        # do not georeference ----
        testthat::expect_warning(
            pr_convert(in_file = testfile, out_folder = out_folder_L1,
                       out_filebase = "testL1_3", out_format = "GTiff",
                       base_georef = FALSE,
                       VNIR = TRUE, selbands_vnir = c(450, 650, 850),
                       SWIR = FALSE, FULL = FALSE, ANGLES = FALSE, PAN = FALSE,
                       CLOUD = FALSE, GLINT = FALSE, LC = FALSE,
                       LATLON = FALSE, ERR_MATRIX = FALSE, apply_errmatrix = FALSE,
                       overwrite = FALSE))

        vnir  <- raster::brick(file.path(out_folder_L1, "testL1_3_HCO_VNIR.tif"))
        # angles retrieved by overwriting with L2 are equal to those of L2
        testthat::expect_equal(dim(vnir), c(1000, 1000, 3),
                               tolerance = 0.001)

        # get HRC and complex ATCOR
        # launch pr_convert to creat VNIR (3 bands), SWIR (1 bands), FULL,
        # ANGLES, LATLON, PAN - save to tiff
        pr_convert(in_file = testfile, out_folder = out_folder_L1,
                   out_filebase = "testL1_1", out_format = "GTiff",
                   VNIR = TRUE, selbands_vnir = c(450),
                   selbands_swir = c(1500), ATCOR = TRUE, ATCOR_wls = c(200,800),
                   SWIR = TRUE, FULL = TRUE, ANGLES = TRUE, PAN = TRUE,
                   CLOUD = TRUE, GLINT = TRUE, LC = TRUE,
                   LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE,
                   overwrite = TRUE, source = "HRC")




    })
