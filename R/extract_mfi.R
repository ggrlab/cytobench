empty_tibble <- tibble::tibble(
    "feature" = NA,
    "negative" = NA,
    "positive" = NA,
    "positive.sd" = NA,
    "negative.sd" = NA
)

#' Extract Single Stain Median Fluorescence Intensity (MFI)
#'
#' This function extracts the median fluorescence intensity (MFI) from single-stain FCS files.
#'
#' @param fcs_dir A character string specifying the directory containing the FCS files. Default is "data-raw/s001".
#' @param regex_singlestain A regular expression pattern to identify single-stain FCS files. Default is "-(CD3-.*)|(none)\\.fcs".
#' @param transform_fun
#' A function to transform the fluorescence values. Default is `function(x) { asinh(x / 1e3) }`.
#' The reported MFIs are calculated as the median of the UNtransformed values. Transformation is
#' only used to cluster the negative and positive populations.
#' @param gating_set_file A gating set file or object to apply to the FCS files. Default is NULL.
#' Usually you want to gate lymphocytes before extracting the MFIs. If a character string is provided,
#' the gating set is loaded from the file. If a (GatingSet) object is provided, it is used directly.
#'
#' @param gate_extract A specific gate to extract from the gating set. Default is NULL.
#' @param multistaining A logical indicating whether to extract multi-stained samples. Default is FALSE.
#' Note that if TRUE, the function will return a dataframe as following:
#'
#' #' \preformatted{
#' # A tibble: 10 x 4
# A tibble: 10 x 7
#'   feature negative positive unstained sample            negative.multi positive.multi
#'   <chr>      <dbl>    <dbl>     <dbl> <chr>                      <dbl>          <dbl>
#' 1 FITC-A      528.    5489.     520.  data-raw/MX.comp~           528.          5489.
#' 2 PE-A        497.    4224.     477.  data-raw/MX.comp~           497.          4224.
#' 3 ECD-A       413.    4260.     411.  data-raw/MX.comp~           413.          4260.
#' 4 PC5.5-A     167.    3306.     145.  data-raw/MX.comp~           167.          3306.
#' 5 PC7-A       134.    3958.     131.  data-raw/MX.comp~           134.          3958.
#' 6 APC-A       120.    1896.     106.  data-raw/MX.comp~           120.          1896.
#' 7 AF700-A     110.    1884.      87.7 data-raw/MX.comp~           110.          1884.
#' 8 AA750-A     290.    3116.     284.  data-raw/MX.comp~           290.          3116.
#' 9 PB-A        422.    4434.     388.  data-raw/MX.comp~           422.          4434.
#' 10 KrO-A       458.    3858.     416.  data-raw/MX.comp~           458.          3858.
#' }
#'
#' If there have been NO single-stained samples (except for the unstained samples),
#' "negative" and "positive" will be filled with the values from the multi-stained samples.
#'
#' @param regex_multistain A regular expression pattern to identify multi-stain FCS files. Default is "(_15-MasterMix)\\.fcs".
#' @param multistain_columns A character vector specifying the relevant columns for multi-staining. Default is c("FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A", "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A").
#'
#' @return A data frame with the extracted MFIs. E.g.:
#' \preformatted{
#' # A tibble: 10 x 4
#'    feature negative positive unstained
#'    <chr>      <dbl>    <dbl>     <dbl>
#'  1 FITC-A      530.   58567.     576.
#'  2 PE-A        526.  149506.     511.
#'  3 ECD-A       453.  107584.     454.
#'  4 PC5.5-A     233.  242867.     166.
#'  5 PC7-A       176.  195307.     138.
#'  6 APC-A       136.  105036.     122.
#'  7 AF700-A     131.   28408.      97.2
#'  8 AA750-A     324.   75247.     296.
#'  9 PB-A        469.   13426.     442.
#' 10 KrO-A       442.    4473.     444.
#' }
#' @examples
#' \dontrun{
#' extracted_mfis <- extract_singlestain_mfi(
#'     "data-raw/s001",
#'     gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
#' )
#' }
#' @export
extract_mfi <- function(fcs_dir = "data-raw/s001",
                        regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
                        transform_fun = function(x) {
                            asinh(x / 1e3)
                        },
                        multistaining = FALSE,
                        regex_multistain = "(_15-MasterMix)\\.fcs",
                        multistain_columns = c(
                            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
                            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
                        ),
                        gating_set_file = NULL,
                        gate_extract = NULL) {
    loaded_fcs <- load_mfi_files(
        fcs_dir = fcs_dir,
        regex_singlestain = regex_singlestain,
        gating_set_file = gating_set_file,
        gate_extract = gate_extract
    )
    if (!multistaining) {
        # Extract the single stainings
        joint_df <- extract_singlestain_mfi_wrapper(loaded_fcs, transform_fun = transform_fun)
    } else {
        # If multistaining is enabled, the following extracts the potential UNSTAINED sample
        joint_df <- extract_singlestain_mfi_wrapper(loaded_fcs, transform_fun = transform_fun, relevant_columns = multistain_columns)

        # Then extract the actually multi-stained sample(s)
        loaded_fcs_multistain <- load_mfi_files(
            fcs_dir = fcs_dir,
            regex_singlestain = regex_multistain,
            gating_set_file = gating_set_file,
            gate_extract = gate_extract
        )
        relevant_mfis_multi <- tryCatch(
            {
                extract_relevant_mfis_multistain(
                    loaded_fcs_multistain,
                    transform_fun = transform_fun,
                    relevant_columns = multistain_columns
                )
            },
            error = function(e) {
                list(tibble::tibble(
                    "sample" = NA,
                    "feature" = NA,
                    "negative" = NA,
                    "positive" = NA,
                    "unstained" = NA,
                    "negative.sd" = NA,
                    "positive.sd" = NA,
                    "unstained.sd" = NA
                ))
            }
        )
        #  Merge single and multi stainings
        joint_df <- dplyr::left_join(
            joint_df,
            relevant_mfis_multi,
            by = "feature",
            suffix = c("", ".multi")
        )
        for (x in c("positive", "negative", "positive.sd", "negative.sd")) {
            if (all(is.na(joint_df[[x]]))) {
                joint_df[[x]] <- joint_df[[paste0(x, ".multi")]]
            }
        }
    }

    # Unname each column
    return(dplyr::mutate(joint_df, dplyr::across(tidyr::everything(), ~ unname(.))))
}

#' Extract Single Stain Median Fluorescence Intensity (MFI)
#' @inheritParams extract_mfi
#' @export
extract_singlestain_mfi <- function(fcs_dir = "data-raw/s001",
                                    regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
                                    transform = function(x) {
                                        asinh(x / 1e3)
                                    },
                                    gating_set_file = NULL,
                                    gate_extract = NULL) {
    .Deprecated("extract_mfi")
    extract_mfi(
        fcs_dir = fcs_dir,
        regex_singlestain = regex_singlestain,
        transform_fun = transform,
        gating_set_file = gating_set_file,
        gate_extract = gate_extract,
        multistaining = FALSE
    )
}

load_mfi_files <- function(fcs_dir = "data-raw/s001",
                           regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
                           gating_set_file = NULL,
                           gate_extract = NULL) {

    # List all files in the specified directory that match the given regex pattern
    dir_files <- list.files(fcs_dir, full.names = TRUE, pattern = regex_singlestain)

    # Load the FCS files into a list of cytosets
    loaded_fcs <- sapply(dir_files, flowWorkspace::load_cytoset_from_fcs, simplify = FALSE)

    # Check if a gating set file is provided
    if (!all(is.null(gating_set_file))) {
        # Load the gating set from file if it's a character string, otherwise use the provided object
        if ("character" %in% class(gating_set_file)) {
            gatingset <- flowWorkspace::load_gs(gating_set_file)
        } else {
            gatingset <- gating_set_file
        }

        # Apply the gating set to each loaded FCS file
        loaded_fcs <- lapply(loaded_fcs, function(ff_x) {
            suppressMessages(
                ff_x_gated <- flowWorkspace::gh_apply_to_cs(gatingset, ff_x, compensation_source = "none")
            )

            # If no specific gate is provided, use the first actual gate (after root)
            if (is.null(gate_extract)) {
                # nolint start
                # [1] "root"
                # [2] "/FITC Control: Control Population"
                # [3] "/FITC Control: Control Population/FITC Control: FITC"
                # [4] "/FITC Control: Control Population/FITC Control: FITC/FITC Control: PE"
                # [5] "/FITC Control: Control Population/FITC Control: FITC/FITC Control: ECD"
                # [6] "/FITC Control: Control Population/FITC Control: FITC/FITC Control: PC5.5"
                # nolint end
                gate_extract <- flowWorkspace::gh_get_pop_paths(ff_x_gated)[2]
                warning("No gate_extract provided, using ", gate_extract)
            }

            # Extract the gated data
            gated_ff <- flowWorkspace::gs_pop_get_data(ff_x_gated, gate_extract) |>
                flowWorkspace::realize_view()
            return(gated_ff)
        })
    }

    return(loaded_fcs)
}

#' Extract Single Stain Median Fluorescence Intensity (MFI)
#' This function extracts the median fluorescence intensity (MFI) from single-stain FCS files.
#'
#' @param loaded_fcs A list of cytosets containing the loaded FCS files.
#' @param transform_fun A function to transform the fluorescence values. Default is `function(x) { asinh(x / 1e3) }`.
#' The reported MFIs are calculated as the median of the UNtransformed values. Transformation is
#' only used to cluster the negative and positive populations.
#' @return A data frame with the extracted MFIs. E.g.:
#' \preformatted{
#' # A tibble: 10 x 4
#'    feature negative positive unstained
#'    <chr>      <dbl>    <dbl>     <dbl>
#'  1 FITC-A      530.   58567.     576.
#'  2 PE-A        526.  149506.     511.
#'  3 ECD-A       453.  107584.     454.
#'  4 PC5.5-A     233.  242867.     166.
#'  5 PC7-A       176.  195307.     138.
#'  6 APC-A       136.  105036.     122.
#'  7 AF700-A     131.   28408.      97.2
#'  8 AA750-A     324.   75247.     296.
#'  9 PB-A        469.   13426.     442.
#' 10 KrO-A       442.    4473.     444.
#' }
extract_singlestain_mfi_wrapper <- function(loaded_fcs,
                                            transform_fun = function(x) {
                                                asinh(x / 1e3)
                                            },
                                            relevant_columns,
                                            ...) {
    # Extract median fluorescence intensities (MFIs) from the loaded FCS files
    relevant_mfis_single <- tryCatch(
        {
            extract_relevant_mfis_singlestain(loaded_fcs, transform_fun = transform_fun, ...)
        },
        error = function(e) {
            list(empty_tibble)
        }
    )
    # Combine single stainings into a single data frame
    single_stainings <- do.call(rbind, relevant_mfis_single[!sapply(relevant_mfis_single, function(x) nrow(x) > 1)])
    if (!missing(relevant_columns)) {
        for (relevant_x in relevant_columns) {
            if (!relevant_x %in% single_stainings$feature) {
                single_stainings <- rbind(single_stainings, tibble::tibble(
                    "feature" = relevant_x,
                    "negative" = NA,
                    "positive" = NA,
                    "negative.sd" = NA,
                    "positive.sd" = NA
                ))
            }
        }
    }
    single_stainings <- single_stainings[!is.na(single_stainings[["feature"]]), ]

    # Extract unstained samples
    relevant_unstained <- relevant_mfis_single[sapply(relevant_mfis_single, function(x) nrow(x) > 1)]

    # Check if there is only one unstained sample
    if (length(relevant_unstained) == 1) {
        names(relevant_unstained) <- "unstained"
    } else {
        warning("More than one unstained sample found, returning each MFI by filename")
    }
    # Join all unstained samples
    all_unstained_joint <- Reduce(dplyr::left_join, relevant_unstained)
    if (!is.null(all_unstained_joint)) {
        # Return the joined data frame of single stainings and unstained samples
        if (any(is.na(single_stainings[["feature"]]))) {
            stop("There are NA values in the feature column of single_stainings, joining with unstained will fail.")
        }
        joint_df <- dplyr::left_join(
            single_stainings,
            all_unstained_joint,
            by = "feature"
        )
    }
    return(joint_df)
}

clustering_seeded_mfi <- function(values, seed, transform_fun, featurename) {
    set.seed(seed)
    clustering <- stats::kmeans(transform_fun(values), centers = 2)
    mfis <- sort(tapply(values, clustering$cluster, median))
    mfis_sd <- sort(tapply(values, clustering$cluster, sd))
    mfis_tib <- tibble::tibble(
        "feature" = featurename,
        "negative" = mfis[1],
        "positive" = mfis[2],
        "negative.sd" = mfis_sd[1],
        "positive.sd" = mfis_sd[2]
    )

    return(mfis_tib)
}
clustering_seeded_mfi_multicolor <- function(values, seed=42, transform_fun, featurename) {
    set.seed(seed)
    clustering <- stats::kmeans(transform_fun(values), centers = 2)
    values_copy <- data.table::data.table(values)
    values_copy[, cluster := clustering$cluster]
    # calculate median per column grouped by cluster with data.table
    medians <- values_copy[, lapply(.SD, median), by = cluster]
    sds <- values_copy[, lapply(.SD, sd), by = cluster]
    index_negative_median <- which.min(apply(medians, 1, sum))
    index_positive_sds <- which.max(apply(sds, 1, sum))

    melted_medians <- data.table::melt(medians, id.vars = "cluster")
    melted_medians[, cluster := ifelse(cluster == index_negative_median, "negative", "positive")]
    cast_medians <- data.table::dcast(melted_medians, variable ~ cluster)

    melted_sds <- data.table::melt(sds, id.vars = "cluster")
    melted_sds[, cluster := ifelse(cluster == index_negative_median, "negative.sd", "positive.sd")]
    cast_sds <- data.table::dcast(melted_sds, variable ~ cluster)

    joint <- dplyr::left_join(
        cast_medians,
        cast_sds,
        by = "variable"
    ) |>
        tibble::as_tibble()
    colnames(joint)[1] <- "feature"
    # # A tibble: 10 × 5
    #    variable     negative positive negative.sd positive.sd
    #    <fct>           <dbl>    <dbl>       <dbl>       <dbl>
    #  1 Blue - 530-A     35.3     332.        68.2       147.
    #  2 Blue - 575-A     25.2     225.        51.9        94.3
    #  3 Blue - 610-A     39.2     449.        94.6       172.
    return(joint)
}

#' Extract Single Stain Median Fluorescence Intensity (MFI)
#' This function extracts the median fluorescence intensity (MFI) from single-stain FCS files.
#' @param loaded_fcs_singlestain A list of cytosets containing the loaded FCS files.
#' @param transform_fun A function to transform the fluorescence values. Default is `function(x) { asinh(x / 1e3) }`.
#' The reported MFIs are calculated as the median of the UNtransformed values. Transformation is only used to cluster the negative and positive populations.
#' @return A data frame with the extracted MFIs. E.g.:
#' \preformatted{
#'   feature negative positive unstained
#'  <chr>      <dbl>    <dbl>     <dbl>
#' 1 FITC-A      530.   58567.     576.
#' 2 PE-A        526.  149506.     511.
#' @export
extract_relevant_mfis_singlestain <- function(loaded_fcs_singlestain,
                                              transform_fun = function(x) {
                                                  x
                                              },
                                              seed = 42) {
    sapply(names(loaded_fcs_singlestain), function(f_x) {
        ff_x <- flowWorkspace::cytoframe_to_flowFrame(loaded_fcs_singlestain[[f_x]][[1]])
        nonempty_channel <- which(flowCore::markernames(ff_x) != "empty")

        # Check if there is more than one non-empty channel
        if (length(nonempty_channel) > 1) {
            warning("More than one non-empty channel in ", f_x)
            return(empty_tibble)
        }

        # If no non-empty channel, it is the unstained sample
        if (length(nonempty_channel) == 0) {
            mfis <- apply(flowCore::exprs(ff_x), 2, median)
            mfis_tib <- tibble::tibble(
                "feature" = names(mfis),
                "unstained" = mfis,
            )
        } else {
            # Cluster the relevant channel into two populations and return the median of both
            values_nonempty <- flowCore::exprs(ff_x)[, names(nonempty_channel)]
            mfis_tib <- clustering_seeded_mfi(
                values = values_nonempty,
                seed = seed,
                transform_fun = transform_fun,
                featurename = names(nonempty_channel)
            )
        }

        return(mfis_tib)
    }, simplify = FALSE)
}

#' Extract Multi-Stain Median Fluorescence Intensity (MFI)
#' This function extracts the median fluorescence intensity (MFI) from multi-stained FCS files.
#' @param loaded_fcs_multistain A list of cytosets containing the loaded FCS files.
#' @param transform_fun A function to transform the fluorescence values. Default is `function(x) { asinh(x / 1e3) }`.
#' The reported MFIs are calculated as the median of the UNtransformed values. Transformation is only used to cluster the negative and positive populations.
#' @param relevant_columns A character vector specifying the relevant columns for multi-staining. Default is c("FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A", "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A").
#' @return A data frame with the extracted MFIs. E.g.:
#' \preformatted{
#'  feature negative positive negative.sd positive.sd
#' <chr>      <dbl>    <dbl>     <dbl>      <dbl>
#' 1 FITC-A      530.   58567.     576.      576.
#' 2 PE-A        526.  149506.     511.      511.
#' @export
extract_relevant_mfis_multistain <- function(loaded_fcs_multistain,
                                             transform_fun = function(x) {
                                                 asinh(x / 1e3)
                                             },
                                             relevant_columns,
                                             seed = 42) {
    if (length(loaded_fcs_multistain) != 1) {
        warning("More than one multistain sample found, returning the MFI of each multistain_column by filename. This is NOT intended!")
    }
    multistain_mfis <- lapply(loaded_fcs_multistain, function(ff_x) {
        values_relevantcols <- flowCore::exprs(flowWorkspace::cytoframe_to_flowFrame(ff_x[, relevant_columns][[1]]))
        clustering_seeded_mfi_multicolor(
            values = values_relevantcols,
            seed = 2,
            transform_fun = function(x) asinh(x / 1e3),
        )
    }) |>
        data.table::rbindlist(idcol = "sample") |>
        tibble::as_tibble()
    return(multistain_mfis)
}
