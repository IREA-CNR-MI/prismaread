context("pr_convert")
skip_on_travis()
skip_on_cran()

test_that("pr_convert works as expected - L1C", {
    skip_on_travis()
    skip_on_cran()
    # VNIR ----
    pr_convert("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1snow/", SWIR = F, overwrite = T, ERR_MATRIX = FALSE, base_georef = F, selbands_vnir = c(750,650,550),
                   CLOUD = T, ATCOR = FALSE, FULL = FALSE, PAN = F, GLINT =F, LC = T, ANGLES = T, LATLON = T)

    # use l2b data for georef
        pr_convert("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                       in_L2_file = "D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/", SWIR = F, overwrite = F, ERR_MATRIX = FALSE, base_georef = T, out_filebase = "l1asl2",
                   CLOUD = T, ATCOR = FALSE, FULL = FALSE, PAN = T, GLINT = T, LC = T, ANGLES = T, LATLON = T)


    # SWIR ----
    pr_convert("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/testL1", SWIR = TRUE, VNIR = FALSE, overwrite = T,base_georef = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # FULL Cube (using already available data) ----
    pr_convert("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/testL1", SWIR = TRUE, VNIR = TRUE, overwrite = FALSE,base_georef = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = TRUE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    pr_convert("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/testL1", SWIR = FALSE, VNIR = FALSE, base_georef = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE, overwrite = T)

})


test_that("pr_convert works as expected - L2B", {
    skip_on_travis()
    skip_on_cran()
    # VNIR ----
    pr_convert("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = FALSE, overwrite = T, LATLON = T, ANGLES = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    pr_convert("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = TRUE, VNIR = FALSE, overwrite = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    pr_convert("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = FALSE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE, overwrite = T)


    # ANGLES and OTHER ----
    pr_convert("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = FALSE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, ANGLES = TRUE,
                   , GLINT = FALSE, overwrite = T)


})

test_that("pr_convert works as expected - L2C", {
    skip_on_travis()
    skip_on_cran()
    # VNIR ----
    pr_convert("D:/prismaread/L2C/PRS_L2C_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2C/testL2C", SWIR = FALSE, overwrite = T, base_georef = TRUE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    pr_convert("D:/prismaread/L2C/PRS_L2C_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2C/testL2C", SWIR = TRUE, VNIR = FALSE, overwrite = T, base_georef = TRUE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    pr_convert("D:/prismaread/L2C/PRS_L2C_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2C/testL2C", SWIR = FALSE, VNIR = FALSE,overwrite = T, base_georef = TRUE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE)

})

test_that("pr_convert works as expected - L2D", {
    skip_on_travis()
    skip_on_cran()
    # VNIR ----
    pr_convert("D:/prismaread/L2D/PRS_L2D_STD_20190616102249_20190616102253_0001.he5",
                   out_format = "ENVI",
                   out_folder = "D:/prismaread/L2D/testL2D/ENVI", SWIR = FALSE,overwrite = T, selbands_vnir = c(750,650,550),
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    pr_convert("D:/prismaread/L2D/PRS_L2D_STD_20190616102249_20190616102253_0001.he5",
                   "D:/prismaread/L2D/testL2D", SWIR = TRUE, VNIR = FALSE,overwrite = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    pr_convert("D:/prismaread/L2D/PRS_L2D_STD_20190616102249_20190616102253_0001.he5",
                   "D:/prismaread/L2D/testL2D_pan", SWIR = FALSE, VNIR = FALSE,overwrite = T,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE)

})
