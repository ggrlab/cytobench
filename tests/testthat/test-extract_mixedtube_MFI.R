test_that("Extract mixed tube MFI", {
    extracted_mfis_singlestain <- extract_mfi(
        "data-raw/MX.compensated-manual.navios.s050/",
        regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
        gating_set_file = "data-raw/MX.compensated-manual.navios.s050/gatingset.s050/MX_p048_s050_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation_navios.flowWorkspace_gatingset"
    )
    # # A tibble: 10 × 4
    #    feature negative positive unstained
    #    <chr>      <dbl>    <dbl>     <dbl>
    #  1 FITC-A      518.   63954.     520.
    #  2 PE-A        550.  164250.     477.
    #  3 ECD-A       467.  124248.     411.
    #  4 PC5.5-A     262.  248894.     145.
    #  5 PC7-A       198.  209206.     131.
    #  6 APC-A       132.  125007.     106.
    #  7 AF700-A     136.   31358.      87.7
    #  8 AA750-A     308.   72735.     284.
    #  9 PB-A        445.   21874.     388.
    # 10 KrO-A       465.    5995.     416.

    extracted_mfis_multistain <- extract_mfi(
        "data-raw/MX.compensated-manual.navios.s050/",
        regex_singlestain = "(none)\\.fcs",
        gating_set_file = "data-raw/MX.compensated-manual.navios.s050/gatingset.s050/MX_p048_s050_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation_navios.flowWorkspace_gatingset",
        multistaining = TRUE,
        regex_multistain = "(_15-MasterMix)\\.fcs",
        multistain_columns = c(
            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
        ),
    )
    #     # A tibble: 10 × 7
    #    feature negative positive unstained sample            negative.multi positive.multi
    #    <chr>      <dbl>    <dbl>     <dbl> <chr>                      <dbl>          <dbl>
    #  1 FITC-A      528.    5489.     520.  data-raw/MX.comp…           528.          5489.
    #  2 PE-A        497.    4224.     477.  data-raw/MX.comp…           497.          4224.
    #  3 ECD-A       413.    4260.     411.  data-raw/MX.comp…           413.          4260.
    #  4 PC5.5-A     167.    3306.     145.  data-raw/MX.comp…           167.          3306.
    #  5 PC7-A       134.    3958.     131.  data-raw/MX.comp…           134.          3958.
    #  6 APC-A       120.    1896.     106.  data-raw/MX.comp…           120.          1896.
    #  7 AF700-A     110.    1884.      87.7 data-raw/MX.comp…           110.          1884.
    #  8 AA750-A     290.    3116.     284.  data-raw/MX.comp…           290.          3116.
    #  9 PB-A        422.    4434.     388.  data-raw/MX.comp…           422.          4434.
    # 10 KrO-A       458.    3858.     416.  data-raw/MX.comp…           458.          3858.

    extracted_mfis_multi_AND_singlestain <- extract_mfi(

        "data-raw/MX.compensated-manual.navios.s050/",
        regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
        gating_set_file = "data-raw/MX.compensated-manual.navios.s050/gatingset.s050/MX_p048_s050_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation_navios.flowWorkspace_gatingset",
        multistaining = TRUE,
        regex_multistain = "(_15-MasterMix)\\.fcs",
        multistain_columns = c(
            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
        ),
    )
    # # A tibble: 10 × 7
    # feature negative positive unstained sample            negative.multi positive.multi
    # <chr>      <dbl>    <dbl>     <dbl> <chr>                      <dbl>          <dbl>
    # 1 FITC-A      518.   63954.     520.  data-raw/MX.comp…           528.          5489.
    # 2 PE-A        550.  164250.     477.  data-raw/MX.comp…           497.          4224.
    # 3 ECD-A       467.  124248.     411.  data-raw/MX.comp…           413.          4260.
    # 4 PC5.5-A     262.  248894.     145.  data-raw/MX.comp…           167.          3306.
    # 5 PC7-A       198.  209206.     131.  data-raw/MX.comp…           134.          3958.
    # 6 APC-A       132.  125007.     106.  data-raw/MX.comp…           120.          1896.
    # 7 AF700-A     136.   31358.      87.7 data-raw/MX.comp…           110.          1884.
    # 8 AA750-A     308.   72735.     284.  data-raw/MX.comp…           290.          3116.
    # 9 PB-A        445.   21874.     388.  data-raw/MX.comp…           422.          4434.
    # 10 KrO-A       465.    5995.     416.  data-raw/MX.comp…           458.          3858.

    testthat::expect_equal(
        extracted_mfis_singlestain,
        extracted_mfis_multi_AND_singlestain[, c("feature", "negative", "positive", "unstained")]
    )

    testthat::expect_true(all(extracted_mfis_multistain[["negative"]] != extracted_mfis_singlestain[["negative"]]))
    testthat::expect_true(all(extracted_mfis_multistain[["positive"]] != extracted_mfis_singlestain[["positive"]]))

    testthat::expect_equal(
        extracted_mfis_multistain[, -c(2:3)],
        extracted_mfis_multi_AND_singlestain[, -c(2:3)]
    )
})
