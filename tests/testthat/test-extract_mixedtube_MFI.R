devtools::load_all()
test_that("Extract mixed tube MFI", {
    extracted_mfis_singlestain <- extract_mfi(
        "data-raw/MX.compensated-manual.navios.s050/",
        regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
        gating_set_file = "data-raw/MX.compensated-manual.navios.s050/gatingset.s050/MX_p048_s050_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation_navios.flowWorkspace_gatingset"
    )
    # # A tibble: 10 × 6
    #    feature negative positive positive.sd negative.sd unstained
    #    <chr>      <dbl>    <dbl>       <dbl>       <dbl>     <dbl>
    #  1 FITC-A      469.   61903.        356.      24023.     520.
    #  2 PE-A        477.  152891.        292.      67811.     477.
    #  3 ECD-A       411.  118859.        395.      48639.     411.
    #  4 PC5.5-A     231.  242080.        748.      90438.     145.
    #  5 PC7-A       156.  200283.        396.      81510.     131.
    #  6 APC-A       117.  121531.        202.      46369.     106.
    #  7 AF700-A     124.   30806.        213.      11375.      87.7
    #  8 AA750-A     287.   71039.        234.      24476.     284.
    #  9 PB-A        413.   21084.        286.       7675.     388.
    # 10 KrO-A       435.    5824.        239.       2034.     416.

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
    # nolint start
    # # A tibble: 10 × 11
    #    feature negative positive negative.sd positive.sd unstained sample                                                                               negative.multi positive.multi positive.sd.multi negative.sd.multi
    #    <chr>      <dbl>    <dbl>       <dbl>       <dbl>     <dbl> <chr>                                                                                         <dbl>          <dbl>             <dbl>             <dbl>
    #  1 FITC-A      528.    5489.       1910.        290.     520.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           528.          5489.              290.             1910.
    #  2 PE-A        497.    4224.       1479.        273.     477.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           497.          4224.              273.             1479.
    #  3 ECD-A       413.    4260.       1685.        255.     411.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           413.          4260.              255.             1685.
    #  4 PC5.5-A     167.    3306.       1549.        223.     145.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           167.          3306.              223.             1549.
    #  5 PC7-A       134.    3958.       1566.        238.     131.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           134.          3958.              238.             1566.
    #  6 APC-A       120.    1896.        691.        184.     106.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           120.          1896.              184.              691.
    #  7 AF700-A     110.    1884.        685.        206.      87.7 data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           110.          1884.              206.              685.
    #  8 AA750-A     290.    3116.       1082.        208.     284.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           290.          3116.              208.             1082.
    #  9 PB-A        422.    4434.       1464.        204.     388.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           422.          4434.              204.             1464.
    # 10 KrO-A       458.    3858.       1238.        222.     416.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_none_I…           458.          3858.              222.             1238.
    # nolint end


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
    # nolint start
    # # A tibble: 10 × 11
    # feature negative positive positive.sd negative.sd unstained sample                                                                         negative.multi positive.multi positive.sd.multi negative.sd.multi
    # <chr>      <dbl>    <dbl>       <dbl>       <dbl>     <dbl> <chr>                                                                                   <dbl>          <dbl>             <dbl>             <dbl>
    # 1 FITC-A      469.   61903.        356.      24023.     520.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           528.          5489.              290.             1910.
    # 2 PE-A        477.  152891.        292.      67811.     477.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           497.          4224.              273.             1479.
    # 3 ECD-A       411.  118859.        395.      48639.     411.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           413.          4260.              255.             1685.
    # 4 PC5.5-A     231.  242080.        748.      90438.     145.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           167.          3306.              223.             1549.
    # 5 PC7-A       156.  200283.        396.      81510.     131.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           134.          3958.              238.             1566.
    # 6 APC-A       117.  121531.        202.      46369.     106.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           120.          1896.              184.              691.
    # 7 AF700-A     124.   30806.        213.      11375.      87.7 data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           110.          1884.              206.              685.
    # 8 AA750-A     287.   71039.        234.      24476.     284.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           290.          3116.              208.             1082.
    # 9 PB-A        413.   21084.        286.       7675.     388.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           422.          4434.              204.             1464.
    # 10 KrO-A       435.    5824.        239.       2034.     416.  data-raw/MX.compensated-manual.navios.s050//MX_p048_s050_comp-manual_ungated_…           458.          3858.              222.             1238.
    # nolint end

    testthat::expect_equal(
        extracted_mfis_singlestain,
        extracted_mfis_multi_AND_singlestain[, c("feature", "negative", "positive", "positive.sd", "negative.sd", "unstained")]
    )

    testthat::expect_true(all(extracted_mfis_multistain[["negative"]] != extracted_mfis_singlestain[["negative"]]))
    testthat::expect_true(all(extracted_mfis_multistain[["negative.sd"]] != extracted_mfis_singlestain[["negative.sd"]]))
    testthat::expect_true(all(extracted_mfis_multistain[["positive"]] != extracted_mfis_singlestain[["positive"]]))
    testthat::expect_true(all(extracted_mfis_multistain[["positive.sd"]] != extracted_mfis_singlestain[["positive.sd"]]))

    testthat::expect_equal(
        extracted_mfis_multistain[, -c(2:5)],
        extracted_mfis_multi_AND_singlestain[, -c(2:5)]
    )
})
