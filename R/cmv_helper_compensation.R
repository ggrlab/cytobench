#' Apply Compensation to CMV Flow Cytometry Data
#'
#' This helper function applies a compensation matrix and autofluorescence correction
#' to CMV flow cytometry data as described in Glehr et al. (TO BE SUBMITTED).
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
#' @keywords relativisation
#' @examples
#' \dontrun{
#'    cytobench::cmv_helper_compensation(
#'        ff_unapplied_compensations,
#'        compensation_type = "manual"
#'    )
#' }
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
