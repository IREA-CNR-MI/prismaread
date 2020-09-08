#' @title pr_check_formulas
#' @description Helper function used to check if formula of a spectral index
#'  can be computed
#' @param indexes either a named list containing indexes names and formulas (
#'  e.g., indexes = list(myindex1 = "R500 / R600",
#'                       myindex2 = "(R800 - R680) / (R800 + R680)")), or a
#'  data.frame containing the columns "Name" and "Formula".
#' @return NULL the function is called for its side effects
#' @details DETAILS
#' @examples
#' \dontrun{
#'  pr_check_formulas(indexes = list(myindex1 = "R500 / R600",
#'                       myindex2 = "(R800 - R680) / (R800 + R680)"))
#' }
#' @rdname pr_check_formulas
#' @export
#' @importFrom stringr str_extract_all
#' @importFrom stats runif
pr_check_formulas <- function(indexes){

    if (!inherits(indexes, "list")) {
        indnames <- as.list(indexes$Name)
        indexes  <- as.list(indexes$Formula)
    } else {
        indnames <- names(indexes)
    }
    if (length(indexes) == 0) {
        stop("`indexes` must be a named list of the kind ",
             "indexes = list(myindex = \"R500\" ,or a data.frame",
             "containing the \"Formula\" and \"Name1\"column. Aborting!")
    }

    getwls <- function(x) {
        stringr::str_extract_all(x, "R[0-9,.]*")[[1]]
    }
    # req_wls <- sort(as.numeric(
    #     unique(unlist((lapply(indexes, FUN = function(x) getwls(x)))))))
    for (ind in seq_along(indexes)) {
        formula <- indexes[[ind]]
        bands   <- getwls(indexes[[ind]])
        indname <- indnames[[ind]]
        for (bb in seq_along(bands)) {
            assign(bands[bb], stats::runif(1))
        }
        try_parse <- try(eval(parse(text = formula)), silent = TRUE)
        if (inherits(try_parse, "try-error")) {
            stop("Check on formula failed for index: ", indname,
                 " Please check your inputs - ",
                 "Aborting!")
        } else {
            message("Checking formula for index: ", indname, " - OK !")
        }
    }
    NULL
}
