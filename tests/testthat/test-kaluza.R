library(testthat)
devtools::load_all()


test_that("Read Kaluza analysis", {
    
    read_analysis <- read_kaluza_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis")
    )
    read_analysis <- read_kaluza_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis"),
        apply_fcs = FALSE
    )
    read_analysis_applied <- read_kaluza_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis"),
        apply_fcs = TRUE,
        verbose = FALSE
    )
})
