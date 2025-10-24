#' Extract Single Stain Median Fluorescence Intensity (MFI)
#'
#' This function extracts the median fluorescence intensity (MFI) from single-stain FCS files.
#'
#' @param fcs_dir A character string specifying the directory containing the FCS files. Default is "data-raw/s001".
#' @param regex_singlestain
#' A regular expression pattern to identify single-stain FCS files. Default is "(-(CD3-.*)|(none))\\.fcs$".
#' @param transform_fun
#' A function to transform the fluorescence values. Default is `function(x) { asinh(x / 1e3) }`.
#' The reported MFIs are calculated as the median of the UNtransformed values. Transformation is
#' only used to cluster the negative and positive populations.
#' @param gating_set_file A gating set file or object to apply to the FCS files. Default is NULL.
#' Usually you want to gate lymphocytes before extracting the MFIs. If a character string is provided,
#' the gating set is loaded from the file. If a (GatingSet) object is provided, it is used directly.
#'
#' @param gate_extract A specific gate to extract from the gating set. Default is NULL.
#'
#' @param regex_multistain
#' A regular expression pattern to identify multi-stain FCS files. Default is "(_15-MasterMix)\\.fcs$".
#' @param multistain_columns
#' A character vector specifying the relevant columns for multi-staining.
#'  Default is c("FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A", "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A").
#' @param ... Additional arguments which are not used here.
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
#'
#'
#' If there are multi-stained samples:
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
#' @export
#' @keywords relativisation
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#'
#' extracted_mfis_singlestain <- extract_mfi(
#'     tmpdir,
#'     regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$"
#' )
#'
#' extracted_mfis_multistain <- extract_mfi(
#'     tmpdir,
#'     regex_singlestain = "(none)\\.fcs$",
#'     multistaining = TRUE,
#'     regex_multistain = "(_15-MasterMix)\\.fcs$",
#'     multistain_columns = c(
#'         "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
#'         "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
#'     ),
#' )
extract_mfi <- function(fcs_dir = "data-raw/s001",
                        regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
                        transform_fun = function(x) {
                            asinh(x / 1e3)
                        },
                        regex_multistain = "(_15-MasterMix)\\.fcs",
                        multistain_columns = c(
                            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
                            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
                        ),
                        gating_set_file = NULL,
                        gate_extract = NULL,
                        ...) {
    joint_df <- tryCatch(
        {
            loaded_fcs <- load_mfi_files(
                fcs_dir = fcs_dir,
                regex_singlestain = regex_singlestain,
                gating_set_file = gating_set_file,
                gate_extract = gate_extract
            )
            # If multistaining is enabled, the following extracts the potential UNSTAINED sample
            extract_singlestain_mfi_wrapper(loaded_fcs, transform_fun = transform_fun, relevant_columns = multistain_columns)
        },
        error = function(e) {
            tibble::tibble(
                "feature" = NA,
                "negative" = NA,
                "positive" = NA,
                "unstained" = NA,
                "negative.sd" = NA,
                "positive.sd" = NA,
                "unstained.sd" = NA
            )
        }
    )

    relevant_mfis_multi <- tryCatch(
        {
            # Then extract the actually multi-stained sample(s)
            loaded_fcs_multistain <- load_mfi_files(
                fcs_dir = fcs_dir,
                regex_singlestain = regex_multistain,
                gating_set_file = gating_set_file,
                gate_extract = gate_extract
            )
            extract_relevant_mfis_multistain(
                loaded_fcs_multistain,
                transform_fun = transform_fun,
                relevant_columns = multistain_columns
            )
        },
        error = function(e) {
            return(tibble::tibble(
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
    if (nrow(relevant_mfis_multi) > 0) {
        #  Merge single and multi stainings
        joint_df <- dplyr::full_join(
            joint_df,
            relevant_mfis_multi,
            by = "feature",
            suffix = c("", ".multi")
        )
    }
    # The following replaces the values of the respective columns
    # by their multi-stained counterparts if they are NA.
    for (x in c("positive", "negative", "positive.sd", "negative.sd")) {
        if (all(is.na(joint_df[[x]]))) {
            joint_df[[x]] <- joint_df[[paste0(x, ".multi")]]
        }
    }

    # Unname each column
    return(dplyr::mutate(joint_df, dplyr::across(tidyr::everything(), ~ unname(.))))
}


#' Extract Single Stain Median Fluorescence Intensity (MFI)
#' @inheritParams extract_mfi
#' @param transform
#' As `transform_fun` in `extract_mfi`.
#' @export
#' @keywords relativisation
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#' loaded_fcs <- load_mfi_files(
#'     fcs_dir = tmpdir,
#'     regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
#'     gating_set_file = NULL,
#'     gate_extract = NULL
#' )
#' extract_relevant_mfis_singlestain(loaded_fcs)
extract_singlestain_mfi <- function(fcs_dir = "data-raw/s001",
                                    regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
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


#' Load and Optionally Gate FCS Files for MFI Extraction
#'
#' This function loads FCS files from a directory that match a specified pattern
#' (typically single-stain or multi-stain), and optionally applies a `GatingSet`
#' to restrict analysis to a specific gated population.
#'
#' @param fcs_dir Character. Directory containing FCS files. Default: `"data-raw/s001"`.
#' @param regex_singlestain Character. Regular expression used to filter FCS filenames.
#'   Only matching files will be loaded. Default: `"(-(CD3-.*)|(none))\\.fcs$"`.
#' @param gating_set_file Optional. A path to a `GatingSet` folder (as a string) or an in-memory `GatingSet` object.
#'   If provided, the gating set will be applied to all loaded FCS samples.
#' @param gate_extract Optional. Name of the population to extract from the gating set.
#'   If `NULL`, the first non-root gate will be used automatically.
#'
#' @return A named list of cytosets or gated cytosets, ready for MFI computation.
#'   Each entry corresponds to one FCS file.
#' @keywords relativisation
#' @keywords internal
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#' loaded_fcs <- load_mfi_files(
#'     fcs_dir = tmpdir,
#'     regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
#'     gating_set_file = NULL,
#'     gate_extract = NULL
#' )
#' @export
load_mfi_files <- function(fcs_dir = "data-raw/s001",
                           regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
                           gating_set_file = NULL,
                           gate_extract = NULL) {
    # List all files in the specified directory that match the given regex pattern
    dir_files <- list.files(fcs_dir, full.names = TRUE)
    dir_files <- dir_files[grepl(regex_singlestain, dir_files)]

    # Load the FCS files into a list of cytosets
    loaded_fcs <- tryCatch(
        {
            sapply(dir_files, flowWorkspace::load_cytoset_from_fcs, simplify = FALSE)
        },
        error = function(e) {
            # In one file (R2_d017_s017_comp-manual_ungated_none_Inf_navios_11-none.fcs),
            # I had the following error:
            #   Error: bad lexical cast: source type value could not be interpreted as target
            # Testing with flowCore::read.FCS gives a warning:
            #   > a <- sapply(dir_files, flowCore::read.FCS, simplify = FALSE)
            # uneven number of tokens: 343
            # The last keyword is dropped.
            # uneven number of tokens: 343
            # The last keyword is dropped.
            # Solution: Rewriting the files into new fcs files.
            warning("Error loading FCS files directly, rewriting them first: ", e$message)
            current_tmpdir <- tempdir()
            dir.create(file.path(current_tmpdir, "fcs_rewrite"), showWarnings = FALSE)
            sapply(dir_files, function(x) {
                ff <- flowCore::read.FCS(x)
                flowCore::write.FCS(ff, filename = file.path(current_tmpdir, "fcs_rewrite", basename(x)))
            })
            res <- sapply(
                file.path(current_tmpdir, "fcs_rewrite", basename(dir_files)),
                flowWorkspace::load_cytoset_from_fcs,
                simplify = FALSE
            )
            unlink(file.path(current_tmpdir, "fcs_rewrite"), recursive = TRUE)
            return(res)
        }
    )

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
#' @param relevant_columns A character vector specifying the relevant columns for MFI extraction.
#' Default is `NA`, which means all columns will be used.
#' @param ... Additional arguments passed to the clustering function.
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
#' @keywords relativisation
#' @keywords internal
#' @export
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#' loaded_fcs <- load_mfi_files(
#'     fcs_dir = tmpdir,
#'     regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
#'     gating_set_file = NULL,
#'     gate_extract = NULL
#' )
#' extract_singlestain_mfi_wrapper(loaded_fcs)
extract_singlestain_mfi_wrapper <- function(loaded_fcs,
                                            transform_fun = function(x) {
                                                asinh(x / 1e3)
                                            },
                                            relevant_columns = NA,
                                            ...) {
    # Extract median fluorescence intensities (MFIs) from the loaded FCS files
    relevant_mfis_single <- tryCatch(
        {
            extract_relevant_mfis_singlestain(loaded_fcs, transform_fun = transform_fun, ...)
        },
        error = function(e) {
            list(empty_tibble())
        }
    )
    # Combine single stainings into a single data frame
    single_stainings <- do.call(rbind, relevant_mfis_single[!sapply(relevant_mfis_single, function(x) nrow(x) > 1)])
    if (!all(is.na(relevant_columns))) {
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
    } else if (length(relevant_unstained) > 1) {
        warning("More than one unstained sample found, returning each MFI by filename")
    }
    # Join all unstained samples
    all_unstained_joint <- Reduce(dplyr::left_join, relevant_unstained)
    if (!is.null(all_unstained_joint)) {
        # Return the joined data frame of single stainings and unstained samples
        if (any(is.na(single_stainings[["feature"]]))) {
            stop("There are NA values in the feature column of single_stainings, joining with unstained will fail.")
        }
        single_stainings <- dplyr::left_join(
            single_stainings,
            all_unstained_joint,
            by = "feature"
        )
    }
    return(single_stainings)
}
#' Clustered MFI Extraction for Single Marker Channel
#'
#' This function performs seeded k-means clustering (with two clusters) on transformed
#' fluorescence values from a single marker channel. It returns the median and standard
#' deviation of the original (untransformed) values per cluster, labeled as negative and positive.
#'
#' @param values Numeric vector. Raw fluorescence values from a single marker channel.
#' @param seed Integer. Random seed to make clustering reproducible.
#' @param transform_fun Function. Transformation to apply before clustering (e.g., `asinh(x / 1e3)`).
#' @param featurename Character. Name of the feature/marker (e.g., `"FITC-A"`), used for labeling output.
#'
#' @return A one-row `tibble` with columns:
#' \describe{
#'   \item{`feature`}{The marker name.}
#'   \item{`negative`, `positive`}{Clustered median intensities (original scale).}
#'   \item{`negative.sd`, `positive.sd`}{Standard deviations within each cluster.}
#' }
#'
#' @keywords relativisation
#' @examples
#' \dontrun{
#' clustering_seeded_mfi(rnorm(100))
#' clustering_seeded_mfi(
#'     rnorm(100),
#'     seed = 123,
#'     transform_fun = function(x) asinh(x / 1e3),
#'     featurename = "Test Marker"
#' )
#' }
clustering_seeded_mfi <- function(values, seed = 42, transform_fun = function(x) {
                                      x
                                  }, featurename = NA) {
    set.seed(seed)

    # Apply transformation before clustering
    clustering <- stats::kmeans(transform_fun(values), centers = 2)

    # Compute median and standard deviation per cluster on original scale
    mfis <- sort(tapply(values, clustering$cluster, stats::median))
    mfis_sd <- sort(tapply(values, clustering$cluster, stats::sd))

    # Format output as a single-row tibble
    mfis_tib <- tibble::tibble(
        "feature" = featurename,
        "negative" = mfis[1],
        "positive" = mfis[2],
        "negative.sd" = mfis_sd[1],
        "positive.sd" = mfis_sd[2]
    )

    return(mfis_tib)
}


#' Clustered MFI Extraction for Multiple Marker Channels
#'
#' This function performs seeded k-means clustering (k = 2) on transformed multi-channel
#' fluorescence values to separate negative and positive cell populations. It computes
#' median, standard deviation, and interquartile range (IQR) of the untransformed
#' values for each marker, returning a tidy `tibble` with one row per marker.
#'
#' The cluster with the lower total median signal across all markers is considered
#' the negative population.
#'
#' @param values A numeric matrix or data.frame of untransformed expression values.
#'   Columns represent marker channels; rows are cells.
#' @param seed Integer. Random seed to ensure reproducibility. Default is `42`.
#' @param transform_fun Function. Transformation to apply prior to clustering
#'   (e.g., `asinh(x / 1e3)`).
#'
#' @return A `tibble` with columns:
#' \describe{
#'   \item{`feature`}{Channel name (from `colnames(values)`)}
#'   \item{`negative`, `positive`}{Clustered median intensities (original scale)}
#'   \item{`negative.sd`, `positive.sd`}{Standard deviations within clusters}
#'   \item{`negative.iqr`, `positive.iqr`}{Interquartile ranges within clusters}
#' }
#'
#' @keywords relativisation
#' @examples
#' \dontrun{
#' df <- data.frame(
#'     marker1 = rnorm(100, mean = 100, sd = 20),
#'     marker2 = rnorm(100, mean = 200, sd = 30)
#' )
#' clustering_seeded_mfi_multicolor(df)
#' clustering_seeded_mfi_multicolor(df, seed = 123, transform_fun = function(x) asinh(x / 1e3))
#' }
clustering_seeded_mfi_multicolor <- function(values,
                                             seed = 42,
                                             transform_fun = function(x) {
                                                 x
                                             }) {
    cluster <- NULL # linting
    set.seed(seed)

    # Cluster cells (rows) using transformed data
    clustering <- stats::kmeans(transform_fun(values), centers = 2)

    # Add cluster assignments to original data
    values_copy <- data.table::data.table(values)
    values_copy[, cluster := clustering$cluster]

    # calculate median per column (marker) grouped by cluster with data.table
    medians <- values_copy[, lapply(.SD, stats::median), by = cluster]

    # Define the "negative" cluster as the one with the lowest sum of medians
    index_negative_median <- which.min(apply(medians, 1, sum))
    clusternumber_negative_median <- medians[["cluster"]][index_negative_median]

    # Format medians
    melted_medians <- data.table::melt(medians, id.vars = "cluster")
    melted_medians[, cluster := ifelse(cluster == clusternumber_negative_median, "negative", "positive")]
    cast_medians <- data.table::dcast(melted_medians, variable ~ cluster)

    # Format standard deviations
    sds <- values_copy[, lapply(.SD, stats::sd), by = cluster]
    melted_sds <- data.table::melt(sds, id.vars = "cluster")
    melted_sds[, cluster := ifelse(cluster == clusternumber_negative_median, "negative.sd", "positive.sd")]
    cast_sds <- data.table::dcast(melted_sds, variable ~ cluster)

    # Format IQRs
    iqrs <- values_copy[, lapply(.SD, stats::IQR), by = cluster]
    melted_iqrs <- data.table::melt(iqrs, id.vars = "cluster")
    melted_iqrs[, cluster := ifelse(cluster == clusternumber_negative_median, "negative.iqr", "positive.iqr")]
    cast_iqrs <- data.table::dcast(melted_iqrs, variable ~ cluster)

    # Join all summaries
    joint <- dplyr::left_join(cast_medians, cast_sds, by = "variable") |>
        dplyr::left_join(cast_iqrs, by = "variable") |>
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


#' Extract MFIs from Single-Stain Cytometry Samples
#'
#' This function extracts the median fluorescence intensity (MFI) for each single-stain sample
#' using seeded k-means clustering (k=2) to separate positive and negative populations. It also identifies
#' unstained samples (i.e., samples with no non-empty channels) and records their median values.
#'
#' @param loaded_fcs_singlestain A named list of cytosets or cytoframes containing single-stain FCS data.
#'   Each element should contain a single sample with one non-empty marker channel, or none (for unstained).
#' The channels must have "empty" as their description if not filled.
#' @param transform_fun Function applied to the data before clustering.
#' Default: `function(x) asinh(x / 1e3)`. Does not affect reported MFIs, only through the clustering.
#' @param seed Integer random seed for reproducible k-means clustering. Default is `42`.
#'
#' @return A named list of `tibble`s, one per input file, each with columns:
#' \describe{
#'   \item{`feature`}{The marker channel name (e.g., `"FITC-A"`).}
#'   \item{`negative`, `positive`}{MFI values per cluster, if applicable.}
#'   \item{`unstained`}{Median intensity if the sample is unstained (no non-empty channels).}
#' }
#'
#' @export
#' @keywords relativisation
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#' loaded_fcs <- load_mfi_files(
#'     fcs_dir = tmpdir,
#'     regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
#'     gating_set_file = NULL,
#'     gate_extract = NULL
#' )
#' extract_singlestain_mfi_wrapper(loaded_fcs)
#' extract_relevant_mfis_singlestain(loaded_fcs)
extract_relevant_mfis_singlestain <- function(loaded_fcs_singlestain,
                                              transform_fun = function(x) {
                                                  x
                                              },
                                              seed = 42) {
    sapply(names(loaded_fcs_singlestain), function(f_x) {
        ff_x <- flowWorkspace::cytoframe_to_flowFrame(loaded_fcs_singlestain[[f_x]][[1]])

        # Identify non-empty marker channels (all empty channels are named "empty")
        nonempty_channel <- which(flowCore::markernames(ff_x) != "empty")

        # Sanity check: more than one non-empty channel is unexpected for single-stain data
        if (length(nonempty_channel) > 1) {
            warning("More than one non-empty channel in ", f_x)
            return(empty_tibble())
        }

        # Handle unstained sample (no active marker channel)
        if (length(nonempty_channel) == 0) {
            mfis <- apply(flowCore::exprs(ff_x), 2, stats::median)
            mfis_tib <- tibble::tibble(
                "feature" = names(mfis),
                "unstained" = mfis
            )
        } else {
            # Cluster the relevant channel into two populations and return the median of both
            # Extract expression values from the single stained channel
            values_nonempty <- flowCore::exprs(ff_x)[, names(nonempty_channel)]

            # Apply seeded k-means to separate negative and positive
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
#'
#' Extracts median fluorescence intensities (MFIs) and associated standard deviations
#' for multiple markers from multi-stained flow cytometry FCS files. This is done by
#' applying k-means clustering (k = 2) on transformed marker intensities to separate
#' negative and positive populations. The MFIs are reported on the original (untransformed) scale.
#'
#' @param loaded_fcs_multistain A list of cytosets containing the loaded multi-stained FCS files.
#' @param transform_fun A transformation function applied prior to clustering. The transformation
#'   helps separate populations. Default is `function(x) asinh(x / 1e3)`. Does not affect reported MFIs.
#' @param relevant_columns A character vector specifying the marker channels to extract
#'   from the multi-stained files. If missing, all columns are used
#' @param seed Integer seed for reproducible k-means clustering. Default is `42`.
#'
#' @return A `tibble` with one row per sample and feature, containing columns:
#' \describe{
#'   \item{`sample`}{File name or sample ID.}
#'   \item{`feature`}{Marker/channel name (e.g., `"FITC-A"`).}
#'   \item{`negative`, `positive`}{Clustered MFI values on original scale.}
#'   \item{`negative.sd`, `positive.sd`}{Standard deviations within clusters.}
#'   \item{`negative.iqr`, `positive.iqr`}{(Optional) Interquartile ranges if computed.}
#' }
#'
#' @export
#' @keywords relativisation
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#' loaded_fcs_multistain <- load_mfi_files(
#'     fcs_dir = tmpdir,
#'     regex_singlestain = "(_15-MasterMix)\\.fcs$",
#'     gating_set_file = NULL,
#'     gate_extract = NULL
#' )
#' extract_relevant_mfis_multistain(
#'     loaded_fcs_multistain
#' )
extract_relevant_mfis_multistain <- function(loaded_fcs_multistain,
                                             transform_fun = function(x) {
                                                 asinh(x / 1e3)
                                             },
                                             relevant_columns = NA,
                                             seed = 42) {
    # Warn if multiple multistained samples were passed (usually only one is expected)
    if (length(loaded_fcs_multistain) > 1) {
        warning("More than one multistain sample found; returning MFIs per sample. This is NOT intended!")
    }

    # Process each multistain sample individually
    multistain_mfis <- lapply(loaded_fcs_multistain, function(ff_x) {
        if (all(is.na(relevant_columns))) {
            relevant_columns <- flowCore::colnames(ff_x[[1]])
        }
        # Extract expression values for relevant columns
        values_relevantcols <- flowCore::exprs(
            flowWorkspace::cytoframe_to_flowFrame(ff_x[, relevant_columns][[1]])
        )

        # Cluster and extract MFIs using k-means
        clustering_seeded_mfi_multicolor(
            values = values_relevantcols,
            seed = seed,
            transform_fun = transform_fun
        )
    }) |>
        data.table::rbindlist(idcol = "sample") |> # Add sample ID column
        tibble::as_tibble()

    return(multistain_mfis)
}
