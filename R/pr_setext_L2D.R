pr_setext_L2D <- function(geo, band) {
    ex   <- matrix(c(geo$xmin - 15,
                     geo$xmin - 15 + dim(band)[2]*30,
                     geo$ymin - 15,
                     geo$ymin - 15 + dim(band)[1]*30),
                   nrow = 2, ncol = 2, byrow = T)
    ex   <- raster::extent(ex)
    band <- raster::setExtent(band, ex, keepres = FALSE)
    return(band)
}
