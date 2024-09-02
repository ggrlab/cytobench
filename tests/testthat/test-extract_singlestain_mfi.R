test_that("Extract single negative MFI", {
    extracted_mfis <- extract_singlestain_mfi(
        "data-raw/s001",
        gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
    )

    # # A tibble: 10 × 4
    #    feature negative positive unstained
    #    <chr>      <dbl>    <dbl>     <dbl>
    #  1 FITC-A      530.   58567.     576.
    #  2 PE-A        526.  149506.     511.
    #  3 ECD-A       453.  107584.     454.
    #  4 PC5.5-A     233.  242867.     166.
    #  5 PC7-A       176.  195307.     138.
    #  6 APC-A       136.  105036.     122.
    #  7 AF700-A     131.   28408.      97.2
    #  8 AA750-A     324.   75247.     296.
    #  9 PB-A        469.   13426.     442.
    # 10 KrO-A       442.    4473.     444.
    testthat::expect_equal(
        extracted_mfis,
        tibble::tibble(
            feature = c("FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A", "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"),
            negative = c(530, 526, 453, 233, 176, 136, 131, 324, 469, 442),
            positive = c(58567, 149506, 107584, 242867, 195307, 105036, 28408, 75247, 13426, 4473),
            unstained = c(576, 511, 454, 166, 138, 122, 97.2, 296, 442, 444)
        ),
        tolerance = .1
    )
})
