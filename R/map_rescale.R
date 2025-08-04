#' Convenience function for rescaling multiple samples
#'
#' Multiple samples can be rescaled at once using this function. It is necessary to
#' supply a "map" of
#'  - which MFI file belongs to which sample. One MFI can belong to multiple samples but also for a sample multiple MFIs can exist.
#' @param sample_dt_list
#'  A named list of data.tables. Each data.table contains the cytometry data for one sample.
#' @param mfi_dt_list
#' A named list of data.tables. Each data.table contains the MFI data for one sample.
#' @param mapping_mfi_to_sample
#' A data.table with two columns: "fileX" and "mfi_file".
#' The "fileX" column contains the sample names and the "mfi_file" column contains the
#' MFI (file) names used to rescale that sample.
#' @param verbose
#' Whether to print progress messages.
#' @param missing_feature
#' If within the MFIs a feature is missing, this feature is rescaled according to the
#' value of this parameter. See \code{\link{rescale}} for details.
#' @param known_missing_features
#'  See \code{\link{rescale}} for details.
#' @param feature_unified_dict
#' See \code{\link{rescale}} for details.
#' @param ...
#' Further arguments to \code{\link{rescale}}.
#'
#' @export
map_rescale <- function(sample_dt_list,
                        mfi_dt_list,
                        mapping_mfi_to_sample,
                        verbose = TRUE,
                        missing_feature = "center_median",
                        known_missing_features = c("FS INT", "FS TOF", "SS INT", "TIME"),
                        feature_unified_dict = feature_unified_dict_default(),
                        sample_identifier_column = "fileX",
                        ...) {
    new_sample_names_unique_check <- c()
    rescaled <- list()

    for (sample_x in names(sample_dt_list)) {
        if (verbose) {
            cat(paste0("Rescaling    \"", sample_x, "\"\n"))
        }
        current_mfi_df <- mapping_mfi_to_sample[mapping_mfi_to_sample[[sample_identifier_column]] == sample_x, ]
        if (nrow(current_mfi_df) == 0) {
            stop(paste0("No MFI file found for \"", sample_x, "\"\n"))
        }

        for (mfi_i in seq_len(nrow(current_mfi_df))) {
            if (nrow(current_mfi_df) == 1) {
                sample_x_newname <- sample_x
            } else {
                sample_x_newname <- paste0(
                    sample_x, "___MFI",
                    ".cd3_", current_mfi_df[mfi_i, ][["lot"]],
                    ".negative_", current_mfi_df[mfi_i, ][["negative"]],
                    ".csv"
                )
            }
            # mfi_mat <- as.matrix(data.table::fread(file = current_mfi_df[mfi_i, ][["mfi_file"]]))
            mfi_mat <- mfi_dt_list[[current_mfi_df[mfi_i, ][["mfi_file"]]]]
            rescaled[[sample_x_newname]] <- rescale(
                sample_to_rescale = sample_dt_list[[sample_x]],
                extracted_mfi = mfi_mat,
                missing_feature = missing_feature,
                feature_unified_dict = feature_unified_dict,
                known_missing_features = known_missing_features,
                ...
            )
            new_sample_names_unique_check <- c(
                new_sample_names_unique_check,
                sample_x_newname
            )
            if (verbose) {
                cat("    rescaled ", sample_x_newname, "\n")
            }
        }
    }
    if (anyDuplicated(new_sample_names_unique_check)) {
        stop("The new sample names are not unique!")
    }
    return(rescaled)
}

#' Creating a mapping for named mfis and read files
#'
#' This is a very specific function to create a mapping between the MFI files and
#' the read files for our usecase.
#'
#' @param named_mfis A named list of MFI files, one element should look like this:
#'
#'
#'  $Align_02.Lot_A.EX
#'                                     FL1 CD3 FITC FL2 CD3 PE FL3 CD3 ECD FL4 CD3 PC5.5 FL5 CD3 PC7 FL6 CD3 APC FL7 CD3 AF700 FL8 CD3 AA750 FL9 CD3 PB FL10 CD3 KrO
#'  MFI POSITIVE POPULATION                32569.91  110339.02   127206.20     532379.81   626938.75   203029.61      61649.41      48720.48   18348.22      7848.84
#'  MFI NEGATIVE POPULATION  UNSTAINED       236.02     213.52      214.93        231.15      190.11      337.84        219.20        234.48     342.80       636.98
#'
#' where the names consists of <sample>.<lot>.<device>
#' @export
#' @return A tibble with the following columns:
#' - mfi_file: The name of the MFI file
#' - negative: The name of the negative control
#' - sample: The name of the sample which should be rescaled with the mfi_file
#' - lot: The lot of the MFI-sample
#' - device: The device used for the sample
#' - lot_verify: The lot of the verification sample
#'
#' @export
create_mapping <- function(named_mfis,
                           sample_names,
                           FUN.sample_names_to_mapping_df = sample_names_to_mapping_df_packageexample,
                           types_single_staining = c("single_CD3", "single_staining")) {
    # usethis::use_pipe() # called that once to enable the pipe operator for the package
    sample_type <- lot <- NULL # linting
    split_mfi_names <- do.call(rbind, strsplit(basename(names(named_mfis)), "\\."))
    mapping_mfi <- tibble::tibble(
        mfi_file = names(named_mfis),
        negative = sub("Negative_", "", basename(dirname(names(named_mfis)))),
        sample = split_mfi_names[, 1],
        lot = split_mfi_names[, 2],
        device = split_mfi_names[, 3]
    )

    samples_mapped <- lapply(sample_names, FUN.sample_names_to_mapping_df)
    mapping_sample <- do.call(rbind, samples_mapped)
    mapped_single <- mapping_sample |>
        dplyr::filter(sample_type %in% types_single_staining) |>
        dplyr::left_join(mapping_mfi, by = c("sample", "lot", "device"), relationship = "many-to-many")
    mapped_single[["lot_verify"]] <- NA

    if (any(is.na(mapped_single[["mfi_file"]]))) {
        print(mapped_single[is.na(mapped_single[["mfi_file"]]), ], n = 10000)
        stop("No MFI file found for the samples shown above in the terminal.")
    }

    mapped_verify <- mapping_sample |>
        dplyr::filter(!sample_type %in% types_single_staining) |>
        dplyr::mutate(lot_verify = lot) |>
        dplyr::select(-lot) |>
        dplyr::left_join(mapping_mfi, by = c("sample", "device"), relationship = "many-to-many")

    mapped_all <- rbind(mapped_single, mapped_verify)
    return(mapped_all)
}

sample_names_to_mapping_df_default <- function(samplename) {
    split_name <- strsplit(samplename, "/")[[1]]
    mapping_sample <- tibble::tibble(
        fileX = samplename,
        sample = sub("[ ].*", "", split_name[, 3]),
        lot = split_name[, 2],
        device = ifelse(grepl("Navios", split_name[, 3], ignore.case = TRUE), "NAVIOS", "EX"),
        sample_type = split_name[, 1]
    )
    return(mapping_sample)
}
#' @export
sample_names_to_mapping_df_packageexample <- function(samplename) {
    split_name <- strsplit(samplename, "/")[[1]]
    mapping_sample <- tibble::tibble(
        fileX = samplename,
        sample = sub("[ ].*", "", basename(samplename)),
        lot = split_name[length(split_name) - 1],
        device = ifelse(
            grepl("Navios", basename(samplename), ignore.case = TRUE),
            "NAVIOS", "EX"
        ),
        sample_type = split_name[length(split_name) - 2]
    )
    return(mapping_sample)
}
