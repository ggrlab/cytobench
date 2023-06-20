# library(testthat)
devtools::load_all()


test_that("Read Kaluza analysis", {
    read_analysis <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis")
    )
    read_analysis <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis"),
        apply_fcs = FALSE
    )
})

test_that("Read and apply Kaluza analysis", {
    read_analysis <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis"),
        apply_fcs = TRUE
    )
})
test_that("Read Kaluza analysis with Kaluza-internal compensation FAILS GENERALLY", {
    # This might be solved by reading the spillover-matrix from the .xml file within the .analysis file
    read_analysis <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "compensation___example_fcs-testing_protocols.analysis"),
        apply_fcs = TRUE,
    )
    # testthat::expect_error(
        kaluza_check_gated(
            read_analysis,
            exported_gates_obj_or_path = file.path("data-raw", "Kaluza_example_fcs", "compensated_example_fcs Ungated Statistics.csv"),
            do_plots = FALSE
        )
    # )
})

test_that("Read Kaluza analysis with gate application", {
    read_analysis_applied <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis"),
        apply_fcs = TRUE,
        verbose = FALSE
    )
    extracted_counts <- tibble::as_tibble(
        flowWorkspace::gs_pop_get_stats(read_analysis_applied$applied_gatings)
    )

    extracted_counts <- extracted_counts %>% dplyr::mutate(
        Gate = stringr::str_replace(pop, ".*/", "")
    )
    extracted_counts$Gate[extracted_counts$Gate == "root"] <- "All"

    kaluza_true_counts <- tibble::as_tibble(data.table::fread(
        file.path("data-raw", "Kaluza_example_fcs", "example_fcs Ungated Statistics.csv"),
        sep = ";", dec = ",",
        skip = 3
    ))

    if (anyDuplicated(kaluza_true_counts$Gate) != 0) {
        stop("Duplicated gate names in example_fcs Ungated Statistics.csv - cannot merge according to them")
    }
    joined_true_mine <- dplyr::full_join(kaluza_true_counts, extracted_counts, by = "Gate")
    different_gates <- joined_true_mine %>%
        dplyr::filter(
            # Number comes from kaluza_true_counts (and is therefore "true")
            # Count comes from extracted_counts (and is therefore "mine")
            Number != count
        ) %>%
        dplyr::mutate(ratio_true_mine = Number / count)
    print(different_gates, n = 10000)
    # gs <- read_analysis_applied$gatingsets[[1]]
    # nodelist <- flowWorkspace::gs_get_pop_paths(gs)
    # str(flowWorkspace::gs_pop_get_gate(gs, "B--")[[1]])
    # str(flowWorkspace::gs_pop_get_gate(gs, "B-+")[[1]])
    # str(flowWorkspace::gs_pop_get_gate(gs, "N+")[[1]])
    different_gates <- different_gates %>%
        dplyr::filter(!Gate %in% c(
            # Exclude the Hinge gates, which I cannot read yet
            "B++", "B--", "B-+", "B+-",
            "H++", "H--", "H-+", "H+-",
            # Exclude the Ellipsoid gates, which I cannot read yet
            "K",
            # Exclude the "O"-gate, which is a linear - still with wrong counts
            # This comes from the fact that within Kaluza, the gate CROSSED the minimum x-axis value
            # and therefore ALL CELLS below were counted being in "O".
            # However, setting the axis properly such that the complete gate is visible
            # returns approximately the same counts here as in Kaluza
            "O"
        ))
    print(different_gates) # The following are the remaining ones. I think they are all correct (enough).
    #     # A tibble: 11 × 9
    # Gate  Number `%Total` `%Gated` Logic       sample      pop      count ratio_true_mine
    # <chr>  <int>    <dbl>    <dbl> <chr>       <chr>       <chr>    <dbl>           <dbl>
    # 1 C       4950    99       99    C           example_fcs /C        4977           0.995
    # 2 F       4955    99.1     99.1  F           example_fcs /F        4959           0.999
    # 3 G--     3483    69.7     69.7  G--         example_fcs /G--      3453           1.01
    # 4 G+-     1134    22.7     22.7  G+-         example_fcs /G+-      1091           1.04
    # 5 I       3973    79.5     79.5  I           example_fcs /I        4009           0.991
    # 6 J       2629    52.6     52.6  J           example_fcs /J        2628           1.00
    # 7 L       4141    82.8     82.8  L           example_fcs /L        4230           0.979
    # 8 N-      1805    36.1     36.1  N-          example_fcs /N-       1806           0.999
    # 9 N+      3195    63.9     63.9  N+          example_fcs /N+       3194           1.00
    # 10 Q--       35     0.7     74.5  Q-- AND P++ example_fcs /P++/Q--    33           1.06
    # 11 Q+-        2     0.04     4.26 Q+- AND P++ example_fcs /P++/Q+-     1           2
})

test_that("Read Kaluza analysis, internal_check", {
    read_analysis <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "analysis___example_fcs-testing_protocols.analysis"),
        apply_fcs = TRUE,
    )

    kaluza_check_gated(
        read_analysis,
        exported_gates_obj_or_path = file.path("data-raw", "Kaluza_example_fcs", "example_fcs Ungated Statistics.csv")
    )
})

test_that("Read Kaluza boolean gates", {
    read_analysis <- kaluza_read_analysis(
        path = file.path("data-raw", "Kaluza_example_fcs", "boolean___example_fcs-testing_protocols.analysis"),
        apply_fcs = FALSE
    )
})
