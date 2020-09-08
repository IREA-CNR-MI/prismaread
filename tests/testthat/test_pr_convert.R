context("Access L2D data")
test_that(
    "Tests on MODIStsp", {
        skip_on_cran()
        skip_on_travis()
        library(piggyback)
        testfile <- file.path(
            usethis::proj_get(),
            "inst/testdata/PRS_L2D_STD_20200524103704_20200524103708_0001.he5")

        # Download and unzip using piggyback if necessary
        if (!file.exists(testfile)){
            pb_download("PRS_L2D_STD_20200524103704_20200524103708_0001.zip",
                        repo = "lbusett/prismaread",
                        dest = file.path(usethis::proj_get(), "inst/testdata"))
            pb_track(glob = "inst/testdata/*.zip, inst/testdata/*.he5")
            unzip(testfile, exdir = dirname(testfile))
            unlink(file.path(
            usethis::proj_get(),
            "inst/testdata/PRS_L2D_STD_20200524103704_20200524103708_0001.zip"))
        }
        out_folder <- file.path(tempdir(), "prismaread/L2D")
        dir.create(out_folder, recursive = TRUE)

        # launch pr_convert to creat VNIR (3 bands), SWIR(bands), FULL,
        # ANGLES, LATLON, PAN
        pr_convert(in_file = testfile, out_folder = out_folder,
                   out_filebase = "testL2D_1", out_format = "GTiff",
                   VNIR = TRUE, selbands_vnir = c(450,550,650,750,850,950),
                   selbands_swir = c(1000,1100,2000,2300),
                   SWIR = TRUE, FULL = TRUE, ANGLES = TRUE,
                   LATLON = TRUE, ERR_MATRIX = TRUE, apply_errmatrix = FALSE)

    })
