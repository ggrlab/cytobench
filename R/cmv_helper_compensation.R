#' CMV compensation helper function
#'
#' This function is a helper function to apply the compensation to the CMV data from
#' Glehr et al. [TO BE SUBMITTED]. It reads the spillover matrix and autofluorescence vector from
#' the FCS file and applies the compensation to the data. Different compensation options are stored
#' within the FCS files and can hereby be selected and applied easily.
#' @param ff_unapplied_compensations
#' A flowFrame object with the unapplied compensations, generated for the CMV data from Glehr et al. [TO BE SUBMITTED]
#' @param compensation_type
#' What compensation to use. Default is "manual". Other options are "auto" and "auto_singlestain".
#' "auto" calculated the spillover matrix and autofluorescence vector  based on the DuraClone
#' Compensation TUBES from Beckman Coulter. "auto_singlestain" calculates the spillover matrix and
#' autofluorescence vector  based on the single stain controls from the CMV data.
#' "manual" is the default option, which uses the manually curated spillover matrix and autofluorescence vector.
#' @param expected_max_proportion_value
#' The expected maximum value for the spillover matrix and autofluorescence vector. Default is 1 which corresponds to 100%. However, technically nothing speaks against spillover values above 100%.
#' Within Kaluza, the spillover values can be a maximum of +/- 500% (+/-5). Thus the default value of 5. This is a safety measure to avoid using the already calculated absolute autofluorescence values.
#' @export
cmv_helper_compensation <- function(ff_unapplied_compensations,
                                    compensation_type = c("manual", "auto", "auto_singlestain"), expected_max_proportion_value = 5) {
    params_spillover_autofluorescence <- read.FCS_custom_spillover(
        ff_unapplied_compensations,
        custom_spillover_keyword = compensation_type[1]
    )
    if (any(
        sapply(params_spillover_autofluorescence, function(x) any(x > expected_max_proportion_value))
    )) {
        stop(
            paste0(
                "Spillover matrix or autofluorescence vector contains values > ",
                expected_max_proportion_value,
                ". This would mean that these are most probably not percentages."
            )
        )
    }
    ff_compensated <- cytoKal::compensation_autofluorescence_kaluza(
        ff_unapplied_compensations,
        spillmat = params_spillover_autofluorescence[["spillover"]],
        autofluorescence_vector = params_spillover_autofluorescence[["autofluorescence_proportional"]],
        add_suffix = "",
        calculate_absolute_values_from_ff = TRUE,
        remove_spillover_keywords = TRUE
    )
    # I checked that the compensation is NOT applied to ff, so I can reuse it
    # head(flowCore::exprs(ff_compensated))
    # head(flowCore::exprs(ff))
    return(ff_compensated)
}


#' Read custom spillover matrix and autofluorescence vector from FCS file
#'
#' Be careful when using it. This is primarily a helper function for the cmv_helper_compensation function and the related data.
#' @param fcs
#' Either a path to an FCS file or a flowFrame object
#' @param custom_spillover_keyword
#' The keyword to use for the custom spillover matrix and autofluorescence vector. Default is "spillover.manual".
#' @export
read.FCS_custom_spillover <- function(fcs, custom_spillover_keyword = "spillover.manual") {
    if (is.character(fcs)) {
        # Reading the stored spillover matrices:
        fs_keywords <- flowCore::read.FCSheader(fcs)[[1]]
        # By reading only the header, the spillover matrix is not transformed into a 
        # matrix. So we have to do it manually.
        fs_keywords["$SPILLOVER"] <- list(
            flowCore:::string_to_spill(fs_keywords[["$SPILLOVER"]])
        )
    } else {
        fs_keywords <- flowCore::keyword(fcs)
    }

    if (!startsWith(custom_spillover_keyword, "spillover.")) {
        custom_spillover_keyword <- paste0("spillover.", custom_spillover_keyword)
    }

    # Kaluza uses the following function to read the spillover matrix:
    # read_spillover <- flowCore:::string_to_spill(fs_keywords[["$SPILLOVER"]])
    current_spillover <- fs_keywords[["$SPILLOVER"]]
    spillovermat <- matrix(
        as.numeric(read.table(text = fs_keywords[[custom_spillover_keyword]])),
        ncol = ncol(current_spillover),
        dimnames = list(
            NULL,
            colnames(current_spillover)
        )
    )
    autofluorescence_proportional <- as.numeric(
        read.table(
            text = fs_keywords[[sub("spillover", "spillover_autofluorescence", custom_spillover_keyword)]]
        )
    )
    names(autofluorescence_proportional) <- colnames(current_spillover)
    return(
        list(
            "spillover" = spillovermat,
            "autofluorescence_proportional" = autofluorescence_proportional
        )
    )
}
