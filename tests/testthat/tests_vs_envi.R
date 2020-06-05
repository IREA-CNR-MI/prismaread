
# L1 ----
in_file <- "D:/prismaread/test/PRS_L1_STD_OFFL_20200524103704_20200524103708_0001.he5"
in_L2_file <- "D:/prismaread/test/PRS_L2B_STD_20200524103704_20200524103708_0001.he5"
out_folder <- "D:/prismaread/test/confronto_envi"

convert_prisma(in_file,
               out_folder, SWIR = F, overwrite = T,
               ERR_MATRIX = FALSE, base_georef = T,
               CLOUD = F, ATCOR = T, FULL = FALSE, PAN = F, GLINT = F, LC = F,
               selbands_vnir = c(500,600,700),
               ANGLES = F, LATLON = F, out_filebase = "432test_l1")

convert_prisma(in_file,
               out_folder, SWIR = F, overwrite = T,
               ERR_MATRIX = FALSE, base_georef = T,
               CLOUD = T, ATCOR = FALSE, FULL = FALSE, PAN = T, GLINT = T, LC = T,
               ANGLES = T, LATLON = T,selbands_vnir = c(500,600,700),
               fill_gaps = FALSE, out_filebase = "432test_l1asl2",
               in_L2_file = in_L2_file)


# L2 ----
in_file <- "D:/prismaread/test/PRS_L2B_STD_20200524103704_20200524103708_0001.he5"
out_folder <- "D:/prismaread/test/confronto_envi"
convert_prisma(in_file,
               out_folder, SWIR = F, overwrite = F,
               ERR_MATRIX = FALSE, base_georef = T,
               CLOUD = F, ATCOR = FALSE, FULL = FALSE, PAN = F, GLINT = F, LC = F,
               ANGLES = F, LATLON = F, selbands_vnir = c(44,34,24),
               out_filebase = "3bands_l2b_nocor")

# L2Dgcp ----
in_file <- "D:/prismaread/test/PRS_L2D_STD_20200523102029_20200523102033_0001.he5"
out_folder <- "D:/prismaread/test/confronto_envi/gcp"
convert_prisma(in_file,
               out_folder, SWIR = F, overwrite = F,
               ERR_MATRIX = FALSE, base_georef = T, VNIR = FALSE,
               CLOUD = F, ATCOR = FALSE, FULL = FALSE, PAN = F, GLINT = F, LC = F,
               ANGLES = T, LATLON = T, selbands_vnir = c(750,650,550),
               out_filebase = "3bands_l2d_nocor")
