# Tests extraction -----

context("Extract statistics")
test_that(
    "Extract statistics", {
        skip_on_cran()
        skip_on_travis()
        #standard, on full image
        in_file <- system.file("/testdata/prismaread_test_HCO_FULL.tif",
                               package = "prismaread")

        if (in_file == "") {
            message("Downloading test data - This may need a long time!")
            piggyback::pb_download(
                "prismaread_test_HCO_FULL.tif",
                repo = "irea-cnr-mi/prismaread",
                dest = file.path(system.file("", package = "prismaread"),
                                 "/testdata"))
        }

        in_vect <- system.file("/extdata/testdata/polys_prisma.gpkg",
                               package = "prismaread")
        test <- pr_extract_spectra(in_file, in_vect, id_field = "field_id")
        testthat::expect_is(test, "data.frame")
        testthat::expect_equal(dim(test), c(231*9, 4))
        testthat::expect_equal(unique(test$ID),
                               c("cloud", "crop_1",  "crop_2",  "forest",
                                 "soil_1",  "soil_2",  "urban",   "urban_2",
                                 "water"))

        #other, on VNIR
        in_file <- system.file("/testdata/prisma_test_HCO_VNIR.tif",
                               package = "prismaread")
        test <- pr_extract_spectra(in_file, in_vect,
                                   stats_format = "wide")
        testthat::expect_is(test, "data.frame")
        testthat::expect_equal(dim(test), c(63, 19))

        test <- pr_extract_spectra(in_file, in_vect,
                                   selstats = c("mean", "stdev", "coeffvar" ,
                                                "min"))
        testthat::expect_equal(dim(test), c(63*9, 6))

        test <- pr_extract_spectra(in_file, in_vect,
                                   selstats = c("mean", "stdev", "coeffvar",
                                                "min"),
                                   quantiles = TRUE ,
                                   percs = c(0.25,0.5,0.75))
        testthat::expect_equal(dim(test), c(63*9, 9))

        test <- pr_extract_spectra(in_file, in_vect, id_field = "field_id",
                                   allpix = TRUE)
        testthat::expect_is(test, "list")

        outfile <- tempfile(fileext = ".RData")
        test <- pr_extract_spectra(in_file, in_vect, id_field = "field_id",
                                   out_file = outfile)
        testthat::expect_true(file.exists(gsub(".RData", "_stats.RData",
                                               outfile)))

        outfile <- tempfile(fileext = ".xlsx")
        test <- pr_extract_spectra(in_file, in_vect, id_field = "field_id",
                                   out_file = outfile)
        testthat::expect_true(file.exists(gsub(".xlsx", "_stats.xlsx",
                                               outfile)))

        outfile <- tempfile(fileext = ".csv")
        test <- pr_extract_spectra(in_file, in_vect, id_field = "field_id",
                                   out_file = outfile)
        testthat::expect_true(file.exists(gsub(".csv", "_stats.csv", outfile)))
    })
