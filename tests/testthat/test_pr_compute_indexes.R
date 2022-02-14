# Tests indexes computation -----

context("Compute spectral indexes")
test_that(
    "Tests on Compute spectral indexes", {
        skip_on_cran()
        skip_on_travis()
        testfile <- file.path(
            system.file("testdata/", package = "prismaread"),
            "PRS_L2D_STD_20200524103704_20200524103708_0001.he5")

        # Download and unzip using piggyback if necessary
        if (!file.exists(testfile)){
            message("Downloading test data - This may need a long time!")
            piggyback::pb_download(
                "PRS_L2D_STD_20200524103704_20200524103708_0001.zip",
                repo = "irea-cnr-mi/prismaread",
                dest = file.path(
                    system.file("", package = "prismaread"), "/testdata"))
            piggyback::pb_track(
                glob = "inst/testdata/*.zip, inst/testdata/*.he5")
            zipfile <- file.path(
                system.file("testdata/", package = "prismaread"),
                "PRS_L2D_STD_20200524103704_20200524103708_0001.zip")
            unzip(zipfile, exdir = dirname(testfile))
            unlink(zipfile)
        }
        out_folder_ind <- file.path(tempdir(), "prismaread/indexes")
        dir.create(out_folder_ind, recursive = TRUE)

        testthat::expect_warning(
            pr_convert(in_file = testfile, out_folder = out_folder_ind,
                       out_filebase = "testL2D_1", out_format = "GTiff",
                       indexes = c("GI", "MSAVI", "kkkk"),
                       overwrite = TRUE, apply_errmatrix = T))
        testthat::expect_true(file.exists(
            file.path(out_folder_ind, "testL2D_1_GI.tif")))

        testthat::expect_true(file.exists(
            file.path(out_folder_ind, "testL2D_1_MSAVI.tif")))

        GI <- raster::brick(file.path(out_folder_ind, "testL2D_1_GI.tif"))
        meangi <- raster::cellStats(GI, "mean", na.rm = TRUE)
        testthat::expect_equal(meangi, 1.046051, tolerance = 0.01)

        MSAVI <- raster::brick(file.path(out_folder_ind, "testL2D_1_MSAVI.tif"))
        meanmsavi <- raster::cellStats(MSAVI, "mean", na.rm = TRUE)
        testthat::expect_equal(meanmsavi, -0.3487982, tolerance = 0.01)

        pr_convert(in_file = testfile, out_folder = out_folder_ind,
                   out_filebase = "testL2D_2", out_format = "GTiff",
                   indexes = c("GI", "MSAVI"),
                   cust_indexes = list(myindex = "R600"),
                   overwrite = TRUE, apply_errmatrix = T)

        myindex <- raster::brick(file.path(out_folder_ind,
                                           "testL2D_2_myindex.tif"))
        meanmyindex <- raster::cellStats(myindex, "mean", na.rm = TRUE)
        testthat::expect_equal(meanmyindex, 0.05779283, tolerance = 0.01)


        testthat::expect_error(
            pr_convert(in_file = testfile, out_folder = out_folder_ind,
                       out_filebase = "testL2D_2", out_format = "GTiff",
                       indexes = c("GI", "MSAVI"),
                       cust_indexes = list(myindex = "R600 - asd"),
                       overwrite = TRUE, apply_errmatrix = T))
    })
