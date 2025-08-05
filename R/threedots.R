#' Programmatic Access to Internal Functions in a Namespace
#'
#' This infix operator provides access to non-exported (internal) functions
#' from a package's namespace, similar to using `pkg:::fun` in base R.
#' Gotten from:
#'  https://stat.ethz.ch/pipermail/r-devel/2013-August/067210.html
#'
#' @param pkg A character string. The name of the package.
#' @param fun A character string. The name of the internal function to access.
#'
#' @return The function object retrieved from the specified package namespace.
#'
#' @details
#' I need this in for the functions defined in the this very script "threedots.R".
#' This is a workaround to access internal functions that are not exported by the package.
#' It allows you to use the syntax `pkg %:::% fun` to retrieve the function `fun` from the package `pkg`.
#'
#' @examples
#' \dontrun{
#' # Access the internal .unitize function from the grid package
#' grid %:::% ".unitize"
#' }
`%:::%` <- function(pkg, fun) {
    # Retrieve an internal (non-exported) function from a package's namespace
    get(fun,
        envir = asNamespace(pkg),
        inherits = FALSE
    )
}

#' string_to_spill
#'
#' This function is a wrapper for the `string_to_spill` function from the `flowCore` package.
#'
#' @param ... Arguments passed to the `string_to_spill` function.
#' @return The result of the `string_to_spill` function.
string_to_spill <- function(...) {
    `%:::%`("flowCore", "string_to_spill")(...)
}
#' spill2txt
#' 
#' This function is a wrapper for the `spill2txt` function from the `flowCore` package.
#' @param ... Arguments passed to the `spill2txt` function.
#' @return The result of the `spill2txt` function.
spill2txt <- function(...) {
    `%:::%`("flowCore", "spill2txt")(...)
}
#' UpdateDerivedValues
#'
#' This function is a wrapper for the `UpdateDerivedValues` function from the `FlowSOM` package.
#' @param ... Arguments passed to the `UpdateDerivedValues` function.
#' @return The result of the `UpdateDerivedValues` function.
UpdateDerivedValues <- function(...) {
    `%:::%`("FlowSOM", "UpdateDerivedValues")(...)
}

#' GeomScattermore
#'
#' This function is a wrapper for the `GeomScattermore` function from the `scattermore` package.
#' @param ... Arguments passed to the `GeomScattermore` function.
#' @return The result of the `GeomScattermore` function.
GeomScattermore <- function(...) {
    `%:::%`("scattermore", "GeomScattermore")
}
