# #' Pseudorelativise a single FCS file
# #'
# #' Computes per-marker MFIs from negative/positive gate pairs and uses them to min-max rescale
# #' the corresponding channels in the raw event matrix.
# #'
# #' This is "pseudorelativisation" because no proper reference sample is used. Instead,
# #' negative and positive populations are defined on the sample itself. This can lead to
# #' artificially removed biological variation if, for example, a marker is globally
# #' upregulated across all cells in a sample.
# #'
# #' @param fcsfile Character(1). Path to an FCS file, a flowFrame or a flowSet
# #' @param file_gatingset Character(1). Path to a FlowJo/flowWorkspace gating set that can be loaded
# #'   via [flowWorkspace::load_gs()].
# #' @param cofactor_namedvec Named numeric vector of cofactors. Currently not used by this function;
# #'   kept for API compatibility with other cytobench workflows.
# #' @param pseudorelativize_gatemap A data.frame/data.table with columns `feature`, `pop_negative`,
# #'   and `pop_positive`. `pop_*` are matched to gate name endings (via [base::basename()]).
# #' For each marker (e.g. FITC_CD16), the corresponding negative and positive gate(s) must be
# #' specified (e.g. "CD16 Negative" and "CD16 Positive" if those are the gate names in the gating set).
# #' If any feature is duplicated, cells from all corresponding gates are pooled to compute MFIs.
# #' @param transform_flowworkspace Optional transform (e.g. a `transformList`) to apply to the
# #'   loaded gating set before gating.
# #' @param new_maxrange Numeric(1) or NULL. If provided, updates all `$PnR` keywords to this value
# #'   after scaling (with safety checks).
# #' @param scalefactor Numeric(1). Multiplicative factor applied after rescaling.
# #'
# #' @return A list with:
# #' - `ff_relativised`: a `flowFrame` with rescaled expression values and updated keywords
# #' - `mfis_relativise`: the per-feature table of negative/positive MFIs used
# #'
# #' @export
# #' @keywords relativisation
# pseudorelativise <- function(fcsfile,
#                              file_gatingset,
#                              pseudorelativize_gatemap,
#                              transform_flowworkspace = NULL,
#                              new_maxrange = 10e6,
#                              scalefactor = 10e3) {
#     stop("In process")
#     if (!is.null(new_maxrange)) {
#         if (!is.numeric(new_maxrange) || length(new_maxrange) != 1 || is.na(new_maxrange) || new_maxrange <= 0) {
#             stop("`new_maxrange` must be a positive numeric scalar or NULL.")
#         }
#     }
#     if (!is.numeric(scalefactor) || length(scalefactor) != 1 || is.na(scalefactor) || scalefactor <= 0) {
#         stop("`scalefactor` must be a positive numeric scalar.")
#     }
#     required_cols <- c("feature", "pop_negative", "pop_positive")
#     if (!all(required_cols %in% colnames(pseudorelativize_gatemap))) {
#         stop(
#             "`pseudorelativize_gatemap` must contain columns: ",
#             paste(required_cols, collapse = ", ")
#         )
#     }

#     # Ensure one row per feature (avoids list-cols after pivot_wider)
#     gatemap_dt <- data.table::as.data.table(pseudorelativize_gatemap)

#     if (is.character(file_gatingset)) {
#         gh <- flowWorkspace::load_gs(file_gatingset)
#     } else {
#         gh <- file_gatingset
#     }
#     if (!is.null(transform_flowworkspace)) {
#         gh <- flowWorkspace::transform(flowWorkspace::gs_clone(gh), transform_flowworkspace)
#     }

#     if (is.character(fcsfile)) {
#         fs <- flowCore::read.flowSet(fcsfile)
#     } else if ("flowFrame" %in% class(fcsfile)) {
#         fs <- flowCore::flowSet(fcsfile)
#     }
#     gated <- cytobench::gate_cells(
#         flowset = fs,
#         gatingset = gh,
#         # gatename is irrelevant here, we just want to extract MFIs which are extracted for all gates anyways
#         gatename = "/",
#         verbose = FALSE,
#         inplace = FALSE
#     )

#     relevant_mfis <- gated$counts |>
#         dplyr::mutate(pop_ending = basename(pop)) |>
#         tidyr::pivot_longer(
#             cols = -dplyr::all_of(c("sample", "pop", "count", "percent_from_parent", "pop_ending")),
#             names_to = "feature",
#             values_to = "mfi"
#         )

#     mfis_relativise <- dplyr::left_join(
#         gatemap_dt |>
#             tidyr::pivot_longer(
#                 cols = c("pop_negative", "pop_positive"),
#                 names_to = "posneg",
#                 values_to = "pop_ending"
#             ),
#         relevant_mfis,
#         by = c("pop_ending", "feature")
#     ) |>
#         tidyr::pivot_wider(
#             id_cols = c("feature"),
#             names_from = "posneg",
#             values_from = c("mfi", "count", "pop")
#         )

#     ff <- flowCore::read.FCS(fcsfile)
#     ff_relativised_exprs <- cytobench::rescale_extracted(
#         sample_to_rescale = flowCore::exprs(ff),
#         extracted_mfi = mfis_relativise,
#         column_negative = "mfi_pop_negative",
#         column_positive = "mfi_pop_positive",
#         scale_column_fun = cytobench::scale_column_minmax
#     )

#     features_present <- intersect(mfis_relativise$feature, colnames(flowCore::exprs(ff)))
#     if (length(features_present) == 0) {
#         stop("None of the `feature`s in `pseudorelativize_gatemap` exist in the FCS expression matrix.")
#     }
#     flowCore::exprs(ff)[, features_present] <-
#         as.matrix(ff_relativised_exprs[, features_present, with = FALSE]) * scalefactor

#     if (!is.null(new_maxrange)) {
#         kws <- flowCore::keyword(ff)
#         maxkw <- kws[grepl("\\\\$P[0-9]+R", names(kws))] |> unlist()
#         maxkw_num <- suppressWarnings(max(as.numeric(maxkw), na.rm = TRUE))
#         if (is.finite(maxkw_num) && maxkw_num > new_maxrange) {
#             stop(
#                 "The previous flowframe max range value (", maxkw_num,
#                 ") exceeds the new_maxrange ", new_maxrange, ". Change new_maxrange."
#             )
#         }

#         exprs_max <- suppressWarnings(max(flowCore::exprs(ff), na.rm = TRUE))
#         if (is.finite(exprs_max) && exprs_max > new_maxrange) {
#             stop(
#                 "The maximum value in the flowframe (", exprs_max,
#                 ") exceeds the new_maxrange ", new_maxrange, ". Change new_maxrange."
#             )
#         }

#         for (x in names(maxkw)) {
#             flowCore::keyword(ff)[[x]] <- sprintf("%d", new_maxrange)
#         }
#     }

#     # Change the $CYT keyword to indicate pseudorelativisation
#     # $CYT usually contains the cytometer name
#     cyt <- flowCore::keyword(ff)[["$CYT"]]
#     if (is.null(cyt) || length(cyt) == 0) {
#         flowCore::keyword(ff)[["$CYT"]] <- "PSEUDORELATIVISED"
#     } else {
#         flowCore::keyword(ff)[["$CYT"]] <- paste0(cyt, "_PSEUDORELATIVISED")
#     }

#     list(
#         ff_relativised = ff,
#         mfis_relativise = mfis_relativise
#     )
# }
