#' @title Compensate a flowframe as Kaluza does
#' @description
#' Comparing this compensation_matrix to Kaluza's compensation matrix shows that the values are the same but transposed.
#'
#' Additionally, Kaluza has an autofluorescence parameter, which seems to be 0.03 for all channels within Kaluza.
#' The documentation says that the value is "0.0309". The 0.03 is only visual and can be seen by changing the number of decimals.
#'
#' In my dissertation (1.4.1.1. Compensation) I describe the compensation workflow. Essentially the formula is:
#'      observed = spillover %*% truth + autofluorescence
#' Then
#'      truth = spillover^-1 %*% (observed - autofluorescence)
#'
#' IN CONTRAST, Kaluza does something different. See
#'      - https://www.beckman.com/search#q=Gallios%20Flow%20Cytometer&t=coveo-tab-techdocs
#'      - B46171AB.pdf: "Protocols/Compensation/Accounting for Autofluorescence"
#' What they do is:
#'  1. Subtract autofluorescence (=324) from the observed data
#'  2. Apply the compensation matrix
#'  3. Add the autofluorescence back
#'
#' Additionally, the autofluorescence value is established by the following (see my email conversation with Beckman Coulter, Ernest Anderson):
#'
#'      Kaluza has an autofluorescence vector available, if you press the button to show it.
#'      The default value that you’ll see is 0.03. If you expand the visible digits you can get
#'      this to 0.0309. How does this relate to the values I shared above? Kaluza (for historical reasons)
#'      chooses to express autofluorescence as a percent. So using the numbers above we can derive the
#'      same value. The standard value, 324 divided by the full range, times 100 gives you the autofluorescence
#'      as a percentage:
#'
#'          A = 324 / 1048575 * 100 = 0.03089907 => rounds to 0.0309
#'
#' @param flowframe A compensated flowframe (flowCore package)
#' @param spillmat
#' Spillover matrix. If a string, it is assumed to be the name of the spillover matrix in the flowframe,
#' extracted by \code{flowCore::spillover(flowframe)[[spillmat]]}.
#' @param autofluorescence_vector
#' A vector of autofluorescence values, one for each channel in the flowframe.
#' Would usually be 0, must be absolute values.
#' @param add_suffix
#' The suffix to add to the autofluorescence vector and the spillover matrix dimension names to match the flowframe.
#' Usually, the autofluorescence vector and the spillover matrix are named without suffix,
#' e.g. FL1, FL2, ..., FLn but the flowframe has suffixes, e.g. FL1-A, FL2-A, ..., FLn-A.
#' @param calculate_absolute_values_from_ff
#' If \code{TRUE}, the autofluorescence vector is multiplied by the maximum value of the flowframe to get the absolute values. Otherwise, the autofluorescence vector is used as is.
#' @param remove_spillover_keywords
#' If \code{TRUE}, the "spillover", "SPILL" and "$SPILLOVER" keywords are removed from the flowframe.
#' Otherwise if the flowframe is later saved to a .fcs file, the (OLD!, not "'spillmat")
#' spillover matrix will be saved to the .fcs file.
#' This is not a problem if the flowframe is only used in R, but a major problem
#' if the flowframe is loaded into Kaluza because then the saved spillover matrix
#' will be used ADDITIONALLY for compensation.
#' @export
#' @examples
# Simulate a basic flowFrame
#' ff <- simulate_ff(columns = paste0("FL", 1:3, "-A"), flowcore = TRUE)
#'
#' # Define a fake spillover matrix (normally extracted from FCS keywords)
#' spill <- matrix(
#'     c(
#'         1, 0.1, 0.05,
#'         0.02, 1, 0.03,
#'         0.01, 0.05, 1
#'     ),
#'     nrow = 3, ncol = 3, byrow = TRUE
#' )
#' # Simulate a spillover matrix which is necessarily positive-definite
#' spill_pd <- base::crossprod(spill, spill) 
#' colnames(spill_pd) <- rownames(spill_pd) <- paste0("FL", 1:3)
#'
#' # Define autofluorescence as percentage of full scale (default: 1048575)
#' autofluorescence_percent <- c(FL1 = 0.0309, FL2 = 0.0309, FL3 = 0.0309)
#'
#' # Apply Kaluza-style compensation
#' ff_kaluza <- compensation_autofluorescence_kaluza(
#'     flowframe = ff,
#'     spillmat = spill_pd,
#'     autofluorescence_vector = autofluorescence_percent,
#'     add_suffix = "-A", # Simulated ff columns are like "FL1-A"
#'     calculate_absolute_values_from_ff = TRUE
#' )
compensation_autofluorescence_kaluza <- function(flowframe,
                                                 spillmat,
                                                 autofluorescence_vector,
                                                 add_suffix = "-A",
                                                 calculate_absolute_values_from_ff = TRUE,
                                                 remove_spillover_keywords = TRUE) {
    colnames(spillmat) <- paste0(colnames(spillmat), add_suffix)
    if (all(is.null(rownames(spillmat)))) {
        rownames(spillmat) <- colnames(spillmat)
    } else {
        rownames(spillmat) <- paste0(rownames(spillmat), add_suffix)
    }
    names(autofluorescence_vector) <- paste0(names(autofluorescence_vector), add_suffix)
    if (calculate_absolute_values_from_ff) {
        autofluorescence_rangemax <- unlist(range(flowframe)["max", names(autofluorescence_vector)])
        autofluorescence_vector_abs <- autofluorescence_rangemax * autofluorescence_vector
    } else {
        autofluorescence_vector_abs <- autofluorescence_vector
    }
    if (sum(autofluorescence_vector_abs) != 0 &&
        sum(autofluorescence_vector_abs) < length(autofluorescence_vector_abs)) {
        warning(
            "The absolute autofluorescence vector seems to be nonzero but very small.\n",
            "The autofluorescence vector is subtracted before compensation and added back after compensation.\n",
            "If you e.g. put the autofluorescence percentage vector from Kaluza, this should actually be not percentages but percent of the maximum value of the Navios device (being 1048575)."
        )
    }

    # 1. Subtract autofluorescence
    ff_mat <- flowCore::exprs(flowframe)[, names(autofluorescence_vector_abs)]
    # ff_mat <- sweep(ff_mat, 2, autofluorescence_vector_abs, "-")
    ff_mat <- sweep_matrix_col(ff_mat, autofluorescence_vector_abs, "-")
    flowCore::exprs(flowframe)[, names(autofluorescence_vector_abs)] <- ff_mat

    # 2. Apply the compensation matrix
    if (!is.matrix(spillmat)) {
        spillmat <- flowCore::spillover(flowframe)[[spillmat]]
    }
    ff_mat_compensated <- ff_mat %*% solve(spillmat)
    # # Usual compensation, with checking:
    # # I do NOT do this, because otherwise I would have to write the AF-subtracted matrix to the flowframe
    # # which is slow!
    # flowframe <- flowCore::compensate(flowframe, spillmat)
    # all(
    #     flowCore::exprs(flowframe)[, names(autofluorescence_vector_abs)] -
    #         ff_mat_compensated[, names(autofluorescence_vector_abs)] == 0
    # )

    if (remove_spillover_keywords) {
        # Remove the spillover and SPILL keywords
        # Otherwise if the flowframe is later saved to a .fcs file, the (OLD!, not "spillmat")
        # spillover matrix will be saved to the .fcs file.
        # This is not a problem if the flowframe is only used in R, but a major problem
        # if the flowframe is loaded into Kaluza because then the saved spillover matrix
        # will be used ADDITIONALLY for compensation.
        for (spillover_keywords in c("spillover", "SPILL", "$SPILLOVER")) {
            flowCore::keyword(flowframe)[[spillover_keywords]] <- NULL
        }
    }

    # 3. Add the autofluorescence back
    ff_mat_compensated <- sweep_matrix_col(ff_mat_compensated, autofluorescence_vector_abs, "+")
    flowCore::exprs(flowframe)[, names(autofluorescence_vector_abs)] <- ff_mat_compensated

    return(flowframe)
}
