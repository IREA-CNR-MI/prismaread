#' @title rastwrite_lines
#' @description Write a raster object by blocks of lines
#' @param rast_in `Raster* object` to be written to disk
#' @param out_file `character` full path of output image
#' @param out_format `character` [\"TIF\" | \"ENVI\"], Default: 'tif'
#' @return the function is called for its side effects
#' @details DETAILS
#' @rdname rastwrite_lines
#' @author Lorenzo Busetto, phD (2017) <lbusett@gmail.com>
#' @importFrom raster nlayers brick raster blockSize writeStart getValues writeValues writeStop

rastwrite_lines <- function(rast_in,
                           out_file,
                           out_format = "tif") {

    if (raster::nlayers(rast_in) > 1) {
        out <- raster::brick(rast_in, values = FALSE)
    } else {
        out <- raster::raster(rast_in)
    }
    bs <-  raster::blockSize(out)

    out <- raster::writeStart(out,
                              filename = out_file,
                              overwrite = TRUE,
                              options = c("COMPRESS=LZW"),
                              datatype = "INT2S")


    for (i in 1:bs$n) {
        v <- raster::getValues(rast_in, row = bs$row[i], nrows = bs$nrows[i] )
        out <- raster::writeValues(out, v, bs$row[i])
    }
    out <- raster::writeStop(out)
    invisible(NULL)
}
