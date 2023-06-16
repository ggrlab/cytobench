test_that("Extract single negative MFI", {
    # devtools::load_all()
    # library(cytobench); library(testthat)
    extracted_mfis <- mfis_from_excel(
        system.file("extdata", "Pre_Arcsinh_MFI_example.xlsx", package = "cytobench"),
        negative_mfi_colname = "MFI NEGATIVE POPULATION  UNSTAINED"
    )
    expect_length(extracted_mfis, 4)

    all_files <- list.files(
        system.file(file.path("extdata", "raw_csv"), package = "cytobench"),
        full.names = TRUE, recursive = TRUE, pattern = "*.csv"
    )
    created_mapping <- create_mapping(
        named_mfis = extracted_mfis,
        sample_names = all_files
    )
    all_files_read <- lapply(all_files, data.table::fread)
    names(all_files_read) <- all_files
    all_files_rescaled <- map_rescale(
        sample_dt_list = all_files_read,
        mfi_dt_list = extracted_mfis,
        mapping_mfi_to_sample = created_mapping
    )

    all_files_read_second <- lapply(all_files, data.table::fread)
    expect_true(all(all_files_read_second[[1]] != all_files_rescaled[[1]]))
    expect_true(all(all_files_read_second[[1]] == all_files_read[[1]]))

    all_files_rescaled_inplace <- map_rescale(
        sample_dt_list = all_files_read,
        mfi_dt_list = extracted_mfis,
        mapping_mfi_to_sample = created_mapping,
        inplace_datatable = TRUE
    )
    expect_true(all(all_files_read_second[[1]] != all_files_read[[1]]))
    expect_true(all(all_files_read[[1]] == all_files_rescaled[[1]]))
})
