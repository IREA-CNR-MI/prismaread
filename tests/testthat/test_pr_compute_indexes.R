# Tests indexes computation -----

context("Compute spectral indexes")
test_that(
    "Tests on L2D", {
        skip_on_cran()
        skip_on_travis()
        in_file <- system.file("inst/testdata/prismaread_test_HCO_FULL.tif",
                               package = "prismaread")
        in_vect <- system.file("inst/extdata/testdata/polys_prisma.gpkg",
                               package = "prismaread")
        test <- pr_extract_spectra(in_file, in_vect, id_field = "field_id")
        testthat::expect_is(test, "data.frame")
        testthat::expect_equal(dim(test), c(231*9, 4))
        testthat::expect_equal(unique(test$ID),
                               c("cloud", "crop_1",  "crop_2",  "forest",
                                 "soil_1",  "soil_2",  "urban",   "urban_2",
                                 "water"))
    })
