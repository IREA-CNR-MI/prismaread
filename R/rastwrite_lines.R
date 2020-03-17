#' @title rastwrite_lines
#' @description Write a raster object by blocks of lines
#' @param rast_in `Raster* object` to be written to disk
#' @param out_file `character` full path of output image
#' @param out_format `character` [\"TIF\" | \"ENVI\"], Default: 'tif'
#' @param proc_lev `character` [\"1\" | \"2D\"], Default: '1'
#' @return the function is called for its side effects
#' @details DETAILS
#' @rdname rastwrite_lines
#' @author Lorenzo Busetto, phD (2017) <lbusett@gmail.com>
#' @importFrom raster nlayers brick raster blockSize writeStart getValues writeValues writeStop

rastwrite_lines <- function(rast_in,
                            out_file,
                            out_format = "tif",
                            proc_lev = "1",
                            scale_min = NULL,
                            scale_max = NULL,
                            join = FALSE) {

    if (raster::nlayers(rast_in) > 1) {
        out <- raster::brick(rast_in, values = FALSE)
    } else {
        out <- raster::raster(rast_in)
    }
    bs <-  raster::blockSize(out)

    if (substring(proc_lev, 1,1) == "1") {
        datatype = "INT2U"
    } else {
        datatype = "FLT4S"
    }
    out <- raster::writeStart(out,
                              filename = out_file,
                              overwrite = TRUE,
                              options = c("COMPRESS=LZW"),
                              datatype = datatype)

    for (i in 1:bs$n) {
        message("Writing Block: ", i, " of: ", bs$n)
        v <- raster::getValues(rast_in, row = bs$row[i], nrows = bs$nrows[i] )
        if (substring(proc_lev, 1, 1) == "2" & !join) {
            v <- scale_min + (v * (scale_max - scale_min)) / 65535
        }
        out <- raster::writeValues(out, v, bs$row[i])
    }
    out <- raster::writeStop(out)
    invisible(NULL)
}
