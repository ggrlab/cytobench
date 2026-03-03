ff_update_names <- function(flowframe, markermap) {
    # match: returns a vector of the positions of (first) matches of its first argument in its second.
    markerorder <- match(flowCore::colnames(flowframe), markermap[["name"]])
    new_cn <- markermap[markerorder, ][["colnames"]]
    new_description <- markermap[markerorder, ][["description_channel.flurochrome.marker.awh"]]
    names(new_description) <- new_cn
    spillmats <- flowCore::spillover(flowframe)

    # Get the old parameternames to map the spillover matrices
    old_parameternames <- flowCore::keyword(flowframe)[paste0("$P", colnames(spillmats[["$SPILLOVER"]]), "N")] |> unlist()
    new_parameternames <- new_cn[match(old_parameternames, markermap[["name"]])]

    # Rename the flowframe
    flowCore::colnames(flowframe) <- new_cn
    flowCore::markernames(flowframe) <- new_description

    # Rename the spillover matrices
    for (spillmat_x in names(spillmats)) {
        if (!is.null(spillmats[[spillmat_x]])) {
            colnames(spillmats[[spillmat_x]]) <- new_parameternames
            rownames(spillmats[[spillmat_x]]) <- NULL
            # test the compensation:
            flowCore::compensate(flowframe[1, ], spillover = spillmats[[spillmat_x]])
            # Replace the spillover matrix in the flowframe keywords
            flowCore::keyword(flowframe)[[spillmat_x]] <- spillmats[[spillmat_x]]
        }
    }

    # Correct autofluorescence keywords
    kw_af <- grep("autofluorescence", names(flowCore::keyword(flowframe)), value = TRUE)
    for (kw_af_single in kw_af) {
        af <- flowCore::keyword(flowframe)[[kw_af_single]]
        names(af) <- new_parameternames
        flowCore::keyword(flowframe)[[kw_af_single]] <- af
    }
    # Cleanup keywords
    flowCore::keyword(flowframe)[is.na(names(flowCore::keyword(flowframe)))] <- NULL
    flowCore::keyword(flowframe)[grepl("flowCore", names(flowCore::keyword(flowframe)))] <- NULL
    return(flowframe)
}
