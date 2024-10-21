test_that("Extract single negative MFI", {
    extracted_mfis <- extract_singlestain_mfi(
        "data-raw/s001",
        gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
    )

    # # A tibble: 10 × 4
    # feature negative positive unstained
    # <chr>      <dbl>    <dbl>     <dbl>
    # 1 FITC-A      542.   60575.     576. 
    # 2 PE-A        551.  162441.     511. 
    # 3 ECD-A       468.  112634.     454. 
    # 4 PC5.5-A     242.  251954.     166. 
    # 5 PC7-A       186.  204156      138. 
    # 6 APC-A       142.  108628.     122. 
    # 7 AF700-A     135.   29050.      97.2
    # 8 AA750-A     333.   77465.     296. 
    # 9 PB-A        478.   14059.     442. 
    # 10 KrO-A       449.    4665.     444.
    testthat::expect_equal(
        extracted_mfis,
        tibble::tibble(
            "feature" = c("FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A", "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"),
            "negative" = c(542, 551, 468, 242, 186, 142, 135, 333, 478, 449),
            "positive" = c(60575, 162441, 112634, 251954, 204156, 108628, 29050, 77465, 14059, 4665),
            "unstained" = c(576, 511, 454, 166, 138, 122, 97.2, 296, 442, 444)
        ),
        tolerance = .1
    )
})
