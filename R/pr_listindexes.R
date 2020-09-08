#' @title pr_listindexes
#' @description helper function used to show the user a DT containing the list
#'  of available spectral indexes
#' @return NULL - The function is called for its side effects
#' @rdname pr_listindexes
#' @export
#' @importFrom DT datatable
pr_listindexes <- function(){
    list <- read.table(system.file("extdata/indexes_list.txt",
                                   package = "prismaread"), sep = "\t",
                       header = T)
    DT::datatable(list)
}
