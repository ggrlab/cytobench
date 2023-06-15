test_that("Extract single negative MFI", {
    extracted_mfis <- mfis_from_excel("inst/extdata/Pre_Arcsinh_MFI_example.xlsx")
    expect_length(extracted_mfis, 4)
})
test_that("Extract multiple negative MFI, unnamed", {
    expect_error({
        extracted_mfis <- mfis_from_excel(
            "inst/extdata/Pre_Arcsinh_MFI_example.xlsx",
            negative_mfi_colname = c(
                "MFI NEGATIVE POPULATION  UNSTAINED",
                "MFI NEGATIVE POPULATION (GG)"
            )
        )
    })
})

test_that("Extract multiple negative MFI, named", {
    extracted_mfis <- mfis_from_excel(
        "inst/extdata/Pre_Arcsinh_MFI_example.xlsx",
        negative_mfi_colname = c(
            "extern" = "MFI NEGATIVE POPULATION  UNSTAINED",
            "intern" = "MFI NEGATIVE POPULATION (GG)"
        )
    )
    expect_length(extracted_mfis, 2)
    expect_length(extracted_mfis[[1]], 4)
    expect_length(extracted_mfis[[2]], 4)
    expect_length(colnames(extracted_mfis[[1]][["Align_02.Lot_A.EX"]]), 10)
    expect_length(rownames(extracted_mfis[[1]][["Align_02.Lot_A.EX"]]), 2)
})
