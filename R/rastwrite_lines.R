rastwrite_lines < function(rast_in,
                           out_file,
                           out_format = "tif") {

    out <- raster::brick(rast_in, values = FALSE)
    bs <-  raster::blockSize(out)
    out_file <- paste(tools::file_path_sans_ext(out_file), "_SWIR")
    out_file <- ifelse(out_format == "tif",
                       paste0(out_file, ".tif"),
                       paste0(out_file, ".envi"))
    out <- raster::writeStart(out, filename = out_file, overwrite = TRUE)

    for (i in 1:bs$n) {
        print(i)
        v <- raster::getValues(rast_swir, row = bs$row[i], nrows = bs$nrows[i] )
        out <- raster::writeValues(out, v, bs$row[i])
    }
    out <- raster::writeStop(out)
}
