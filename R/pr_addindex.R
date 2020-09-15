#' @title add a spectral index to the list of available ones
#' @description Function used to add a spectral index to the list of available
#' ones
#' @param Name `character` Name of the new index (e.g., "myindex")
#' @param Formula `character` Formula of the index. In the formula, bands
#'  are to be referred to the prefix "R", followed by the required wavelength
#'  /e.g., "R600 / R720"
#' @param Description `character`, Optional description of the index (e.g.,
#' "My Index for Nitrogen"), Default: NULL
#' @param Sensitivity `character`, Optional specification of index sensitivity
#'   (e.g., "Nitrogen concentration"),  Default: NULL
#' @param Reference `character`, Optional specification of index literature
#'  reference (e.g., "me et al., (2020)), Default: NULL
#' @return NULL - The function is called for its side effects
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#' pr_addindex(Name = "myindex", Formula = "R600 / R700",
#'            Description = "My custom Index",
#'            Reference = "Me (2020)")
#'  }
#' }
#' @rdname pr_addindex

pr_addindex <- function(Name,
                        Formula,
                        Description = NA,
                        Sensitivity = NA,
                        Reference   = NA) {
    curr_indexes <- read.table(system.file("extdata/indexes_list.txt",
                                           package = "prismaread"), sep = "\t",
                               header = TRUE)
    newindex <- data.frame("Name"        = Name,
                           "Description" = Description,
                           "Formula"     = Formula,
                           "Sensitivity" = Sensitivity,
                           "Reference"   = Reference,
                           stringsAsFactors = FALSE)
    if (newindex$Name %in% curr_indexes$Name) {
        stop("A index with the same list alrasy exists in the list. Aborting!")
    } else {
        newindexes <- rbind(curr_indexes, newindex)
        write.table(newindexes,
                    file = system.file("extdata/indexes_list.txt",
                                       package = "prismaread"), sep = "\t",
                    col.names = TRUE)
    }
    invisible(NULL)
}
