test_that("Extract single and mixed tube MFI", {
    extracted_mfis <- extract_mfi(
        "testdata_nogit/2025-09-17_s017NoMFIextract/s017",
        regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
        transform = function(x) {
            asinh(x / 5e3)
        },
        # If no gating_set_file is given, use all cells.
        # In this case though, gating happened previously.
        # gating_set_file = gatingset,
        # gate_extract = gate_extract,
        filepattern = paste0("((_",
            # 12: panel
            # 13: singlestaining+unstained joint
            # 14: singlestaining+unstained joint, subsampled to 100k cells  -> This is AGAIN subsampled to 10k, so redundant to 13 (after subsampling to 10k!)
            sprintf("%02d", c(1:12, 14)),
            "-).*\\.fcs$)",
            collapse = "|"
        ),
        column_positive = "positive",
        column_negative = "negative",
    )
})
