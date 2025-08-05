test_that("Extract single and mixed tube MFI", {
    set.seed(42)
    fs_ss <- simulate_cd3()
    tmpdir <- withr::local_tempdir()
    flowCore::write.flowSet(fs_ss, tmpdir)

    extracted_mfis_singlestain <- extract_mfi(
        tmpdir,
        regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$"
    )

    extracted_mfis_multistain <- extract_mfi(
        tmpdir,
        regex_singlestain = "(none)\\.fcs$",
        multistaining = TRUE,
        regex_multistain = "(_15-MasterMix)\\.fcs$",
        multistain_columns = c(
            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
        ),
    )
    # nolint start
    # Browse[1]> extracted_mfis_multistain
    # # A tibble: 10 x 13
    #    feature negative positive negative.sd positive.sd unstained sample
    #    <chr>      <dbl>    <dbl>       <dbl>       <dbl>     <dbl> <chr>
    #  1 FITC-A      10.1 102679.        0.811      10844.   0.171   /tmp/Rtmp3hskjP/~
    #  2 PE-A     80449.   87492.    49878.         52280.  -0.103   /tmp/Rtmp3hskjP/~
    #  3 ECD-A    36915.   84650.    50891.         50487.   0.135   /tmp/Rtmp3hskjP/~
    #  4 PC5.5-A     11.6  92260.    51044.         53664.  -0.0573  /tmp/Rtmp3hskjP/~
    #  5 PC7-A    87958.      10.8   47900.         51403.  -0.00692 /tmp/Rtmp3hskjP/~
    #  6 APC-A       11.3  90162.    49777.         49448.  -0.156   /tmp/Rtmp3hskjP/~
    #  7 AF700-A  93295.      11.2   51885.         49657.  -0.216   /tmp/Rtmp3hskjP/~
    #  8 AA750-A     11.0  42967.    51182.         52945.  -0.0440  /tmp/Rtmp3hskjP/~
    #  9 PB-A     86776.      11.1   52116.         49031.   0.00720 /tmp/Rtmp3hskjP/~
    # 10 KrO-A       12.0     11.9   51716.         51720.   0.0390  /tmp/Rtmp3hskjP/~
    # # i 6 more variables: negative.multi <dbl>, positive.multi <dbl>,
    # #   negative.sd.multi <dbl>, positive.sd.multi <dbl>, negative.iqr <dbl>,
    # #   positive.iqr <dbl>
    # nolint end

    extracted_mfis_multi_AND_singlestain <- extract_mfi(
        tmpdir,
        regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
        multistaining = TRUE,
        regex_multistain = "(_15-MasterMix)\\.fcs$",
        multistain_columns = c(
            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
        )
    )
    # Browse[1]> extracted_mfis_multi_AND_singlestain
    # # A tibble: 10 x 13
    # feature negative positive negative.sd positive.sd unstained sample
    # <chr>      <dbl>    <dbl>       <dbl>       <dbl>     <dbl> <chr>
    # 1 FITC-A      9.97  100470.       1.03        8353.   0.171   /tmp/Rtmp3hskjP/~
    # 2 PE-A       10.2   100541.       1.11       10224.  -0.103   /tmp/Rtmp3hskjP/~
    # 3 ECD-A       9.92  100467.       1.10        9321.   0.135   /tmp/Rtmp3hskjP/~
    # 4 PC5.5-A     9.93   98747.       0.967       9995.  -0.0573  /tmp/Rtmp3hskjP/~
    # 5 PC7-A       9.73   98806.       0.932       8938.  -0.00692 /tmp/Rtmp3hskjP/~
    # 6 APC-A      10.3    98064.       0.903       9981.  -0.156   /tmp/Rtmp3hskjP/~
    # 7 AF700-A     9.92  100866.       1.00       11754.  -0.216   /tmp/Rtmp3hskjP/~
    # 8 AA750-A     9.81   99310.       1.08        8817.  -0.0440  /tmp/Rtmp3hskjP/~
    # 9 PB-A       10.0    99647.       0.933       8512.   0.00720 /tmp/Rtmp3hskjP/~
    # 10 KrO-A       9.75   98631.       0.985      10609.   0.0390  /tmp/Rtmp3hskjP/~
    # # i 6 more variables: negative.multi <dbl>, positive.multi <dbl>,
    # #   negative.sd.multi <dbl>, positive.sd.multi <dbl>, negative.iqr <dbl>,
    # #   positive.iqr <dbl>

    testthat::expect_equal(
        extracted_mfis_singlestain[, c("feature", "negative", "positive", "positive.sd", "negative.sd", "unstained")],
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


test_that("Extract mixed tube MFI, new function", {
    set.seed(42)
    fs_ss <- simulate_cd3()
    tmpdir <- withr::local_tempdir()
    flowCore::write.flowSet(fs_ss, tmpdir)


    extracted_mfis_multistain <- extract_mfi(
        tmpdir,
        regex_singlestain = "(none)\\.fcs$",
        multistaining = TRUE,
        regex_multistain = "(_15-MasterMix)\\.fcs$",
        multistain_columns = c(
            "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
            "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
        ),
        # This transform function
        #   1) is ONLY used for the multistain sample
        #   2) is only relevant for the CLUSTERING; the MFIs are reported without transformations!
        transform_fun = function(x) {
            asinh(x)
        },
    )
    testthat::expect_true(TRUE) # just a run-through test
})
