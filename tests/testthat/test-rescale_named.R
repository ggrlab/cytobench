test_that("Testing rescale_extracted", {
    extracted_mfis <- extract_singlestain_mfi(
        "data-raw/s001",
        gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
    )
    loaded_fcs <- flowCore::read.FCS("data-raw/s001/CT_p001_s001_comp-manual_ungated_none_100k_navios_14-01..11-SingleNone.fcs")
    rescaled_sample <- rescale_extracted(
        sample_to_rescale = flowCore::exprs(loaded_fcs),
        extracted_mfi = extracted_mfis
    )
    testthat::expect(TRUE)
})
