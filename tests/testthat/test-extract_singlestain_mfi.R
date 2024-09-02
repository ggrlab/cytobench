devtools::load_all()

test_that("Extract single negative MFI", {
    extracted_mfis <- extract_singlestain_mfi(
        "data-raw/s001",
        gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
    )
})
