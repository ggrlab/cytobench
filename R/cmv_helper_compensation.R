#' Apply Compensation to CMV Flow Cytometry Data
#'
#' This helper function applies a compensation matrix and autofluorescence correction
#' to CMV flow cytometry data as described in Glehr et al. [TO BE SUBMITTED].
#' It retrieves compensation parameters from the FCS file metadata using a custom keyword.
#' Multiple compensation types can be selected depending on how the compensation
#' was estimated (manual, auto, or auto_singlestain).
#'
#' @param ff_unapplied_compensations
#' A `flowFrame` object containing uncompensated data and
#'   embedded compensation parameters (spillover matrix and autofluorescence vector) in the metadata.
#' @param compensation_type
#' Character. Specifies the type of compensation to apply.
#'   Options are:
#'   \describe{
#'     \item{"manual"}{Use a manually curated spillover matrix and autofluorescence vector (default).}
#'     \item{"auto"}{Use values computed from DuraClone compensation tubes (Beckman Coulter).}
#'     \item{"auto_singlestain"}{Use values derived from single stain controls within the CMV dataset.}
#'   }
#' @param expected_max_proportion_value
#' Numeric. Upper bound for expected spillover and autofluorescence values.
#'   This prevents interpreting absolute values as proportions. Defaults to 5 (corresponding to 500%).
#' @param original_spillover_keyword
#' Character. The keyword in the FCS file that stores the original spillover matrix.
#'   Defaults to `"$SPILLOVER"`.
#'
#' @return A `flowFrame` with compensation and autofluorescence correction applied.
#' @export
cmv_helper_compensation <- function(ff_unapplied_compensations,
                                    compensation_type = c(
                                        "manual", "auto",
                                        "auto_singlestain"
                                    ),
                                    expected_max_proportion_value = 5,
                                    original_spillover_keyword = "$SPILLOVER") {
    # Read the spillover matrix and autofluorescence vector from FCS file using selected compensation type
    params_spillover_autofluorescence <- read.FCS_custom_spillover(
        fcs = ff_unapplied_compensations,
        custom_spillover_keyword = compensation_type[1],
        original_spillover_keyword = original_spillover_keyword
    )

    # Sanity check: ensure values are in expected proportional range (e.g., 0–1 or 0–5)
    if (any(
        sapply(params_spillover_autofluorescence, function(x) any(x > expected_max_proportion_value))
    )) {
        stop(
            paste0(
                "Spillover matrix or autofluorescence vector contains values > ",
                expected_max_proportion_value,
                ". This likely indicates that values are absolute, not proportions."
            )
        )
    }

    # Apply spillover and autofluorescence compensation
    ff_compensated <- cytoKal::compensation_autofluorescence_kaluza(
        ff_unapplied_compensations,
        spillmat = params_spillover_autofluorescence[["spillover"]],
        autofluorescence_vector = params_spillover_autofluorescence[["autofluorescence_proportional"]],
        add_suffix = "", # Keep channel names unchanged
        calculate_absolute_values_from_ff = TRUE,
        remove_spillover_keywords = TRUE
    )

    # Note: Original ff is not compensated — safe to reuse
    return(ff_compensated)
}


#' Read custom spillover matrix and autofluorescence vector from FCS file
#'
#' Be careful when using it. This is primarily a helper function for the cmv_helper_compensation function and the related data.
#' @param fcs
#' Either a path to an FCS file or a flowFrame object
#' @param custom_spillover_keyword
#' The keyword to use for the custom spillover matrix and autofluorescence vector. Default is "spillover.manual".
#' @param original_spillover_keyword
#' The keyword to use for the original spillover matrix. Default is "$SPILLOVER".
#' @return
#' A list with the spillover matrix and the autofluorescence vector
#' @export
read.FCS_custom_spillover <- function(fcs, custom_spillover_keyword = "spillover.manual", original_spillover_keyword = "$SPILLOVER") {
    if (is.character(fcs)) {
        # Reading the stored spillover matrices:
        fs_keywords <- flowCore::read.FCSheader(fcs)[[1]]
        # By reading only the header, the spillover matrix is not transformed into a
        # matrix. So we have to do it manually.
        fs_keywords[original_spillover_keyword] <- list(
            flowCore:::string_to_spill(fs_keywords[[original_spillover_keyword]])
        )
    } else {
        fs_keywords <- flowCore::keyword(fcs)
    }

    if (!startsWith(custom_spillover_keyword, "spillover.")) {
        custom_spillover_keyword <- paste0("spillover.", custom_spillover_keyword)
    }

    # Kaluza uses the following function to read the spillover matrix:
    # read_spillover <- flowCore:::string_to_spill(fs_keywords[[original_spillover_keyword]])

    spillovermat <- flowCore:::string_to_spill(fs_keywords[[custom_spillover_keyword]])
    if (custom_spillover_keyword == "spillover.original") {
        autofluorescence_proportional <- spillovermat[1, ]
        autofluorescence_proportional[TRUE] <- 0
    } else {
        autofluorescence_proportional <- as.numeric(
            read.table(text = fs_keywords[[sub(
                "spillover",
                "spillover_autofluorescence", custom_spillover_keyword
            )]])
        )
        names(autofluorescence_proportional) <- colnames(spillovermat)
    }
    current_spillover <- fs_keywords[[original_spillover_keyword]]
    # Order the spillover matrix according to the current spillover matrix
    spillover_reordered <- spillovermat
    rownames(spillover_reordered) <- colnames(spillover_reordered)
    spillover_reordered <- spillover_reordered[colnames(current_spillover), colnames(current_spillover)]
    return(
        list(
            spillover = spillover_reordered,
            autofluorescence_proportional = autofluorescence_proportional[colnames(current_spillover)]
        )
    )
}
