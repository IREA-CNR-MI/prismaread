context("convert_prisma")


test_that("convert_prisma works as expected - L1C", {

    # VNIR ----
    convert_prisma("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/testL1", SWIR = FALSE, overwrite = T, ERR_MATRIX = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    convert_prisma("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/testL1", SWIR = TRUE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    convert_prisma("D:/prismaread/L1/PRS_L1_STD_OFFL_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L1/testL1", SWIR = FALSE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE)

})


test_that("convert_prisma works as expected - L2B", {

    # VNIR ----
    convert_prisma("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    convert_prisma("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = TRUE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    convert_prisma("D:/prismaread/L2B/PRS_L2B_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2B/testL2B", SWIR = FALSE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE)


})

test_that("convert_prisma works as expected - L2C", {

    # VNIR ----
    convert_prisma("D:/prismaread/L2C/PRS_L2C_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2C/testL2C", SWIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    convert_prisma("D:/prismaread/L2C/PRS_L2C_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2C/testL2C", SWIR = TRUE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    convert_prisma("D:/prismaread/L2C/PRS_L2C_STD_20200215103028_20200215103033_0001.he5",
                   "D:/prismaread/L2C/testL2C", SWIR = FALSE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE)


})

test_that("convert_prisma works as expected - L2D", {

    # VNIR ----
    convert_prisma("D:/prismaread/L2D/PRS_L2D_STD_20190616102249_20190616102253_0001.he5",
                   "D:/prismaread/L2D/testL2D", SWIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # SWIR ----
    convert_prisma("D:/prismaread/L2D/PRS_L2D_STD_20190616102249_20190616102253_0001.he5",
                   "D:/prismaread/L2D/testL2D", SWIR = TRUE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = FALSE, GLINT = FALSE)

    # PAN ----
    convert_prisma("D:/prismaread/L2D/PRS_L2D_STD_20190616102249_20190616102253_0001.he5",
                   "D:/prismaread/L2D/testL2D", SWIR = FALSE, VNIR = FALSE,
                   CLOUD = FALSE, ATCOR = FALSE, FULL = FALSE, PAN = TRUE
                   , GLINT = FALSE)

})
