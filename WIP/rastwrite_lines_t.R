#' @title rastwrite_lines
#' @description Write a raster object by blocks of lines
#' @param rast_in `Raster* object` to be written to disk
#' @param out_file `character` full path of output image
#' @param out_format `character` [\"TIF\" | \"ENVI\"], Default: 'tif'
#' @param proc_lev `character` [\"1\" | \"2D\"], Default: '1'
#' @param scale_min `numeric` coefficients use to compute values from DN on L2 products
#' @param scale_max `numeric`  coefficients use to compute values from DN on L2 products
#' @param join `logical` flag used to indicate if we are saving the "joined" VNIR+SWIR cube
#' @return the function is called for its side effects
#' @details DETAILS
#' @rdname rastwrite_lines
#' @author Lorenzo Busetto, phD (2017)
#' @importFrom raster nlayers brick raster blockSize writeStart getValues writeValues writeStop

rastwrite_lines_t <- function(rast_in,
                              out_file,
                              out_format = "tif",
                              proc_lev = "1",
                              scale_min = NULL,
                              scale_max = NULL,
                              join = FALSE) {

    rast_in <- terra::rast(rast_in)
    out <- terra::rast(rast_in)


    if (proc_lev == "ERR") {
        datatype = "INT1U"
    }

    if (substring(proc_lev, 1,1) == "1") {
        datatype = "FLT4S"
    } else {
        datatype = "FLT4S"
    }

    outbs <- terra::writeStart(out,
                             filename = out_file,
                             overwrite = TRUE,
                             wot = list(gdal = c("COMPRESS=LZW"),
                                        datatype = datatype)
    )

    for (i in 1:outbs$n) {
        message("Writing Block: ", i, " of: ", outbs$n)
        v <- terra::readValues(rast_in, row = outbs$row[i], nrows = outbs$nrows[i] )
        if (substring(proc_lev, 1, 1) == "2" & !join) {
            v <- scale_min + (v * (scale_max - scale_min)) / 65535
        }
        terra::writeValues(out, v, outbs$row[i], nrows = outbs$nrows[i])
    }
    out <- terra::writeStop(out)

    invisible(NULL)
}
