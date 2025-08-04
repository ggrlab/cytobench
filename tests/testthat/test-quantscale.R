test_that("Testing extract_singlestain_quantilefun", {
    test_datadir <- file.path(testthat::test_path(), "testdata", "gated_lymphocytes")
    relevant_cols <- c(
        "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
        "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
    )
    loaded_fcs <- cytobench:::load_mfi_files(
        fcs_dir = file.path(test_datadir, "navios"),
        regex_singlestain = "-(CD3-.*)|(none)\\.fcs"
    )
    qfuns <- extract_singlestain_quantilefun(
        loaded_fcs,
        transform_fun = function(x) {
            asinh(x / 5e3)
        }, relevant_columns = relevant_cols
    )

    # The following shows the generated quantile functions (x=quantile, y=respective MFI)
    # and we observe that the "Quantile-MFI" actually deviates considerably from the
    # previously used MFI.
    # Our base relativisation approach would have "taken" the median fluorescence intensity
    # of the whole positive population, (the 50% quantile).
    quantiles <- seq(0, 1, length.out = 103)
    plotting_df <- lapply(names(unlist(qfuns)), function(x) {
        data.table::data.table(
            marker = x,
            quantiles = quantiles,
            approximation = unlist(qfuns)[[x]](quantiles)
        ) |> tidyr::pivot_longer(
            cols = -c("quantiles", "marker"),
            names_to = "cluster",
            values_to = "value"
        )
    }) |>
        data.table::rbindlist() |>
        dplyr::mutate(
            markername = stringr::str_remove(marker, ".*\\.fcs\\.")
        )
    pdf("removeme.pdf", width = 12)
    print(
        ggplot2::ggplot(plotting_df |> dplyr::filter(
            !grepl("TIME", marker),
            !grepl("SC", marker),
        ), ggplot2::aes(x = quantiles, y = value, col = markername, group = marker)) +
            ggplot2::geom_line()  +
            ggplot2::scale_y_log10()
    )
    dev.off()
})

# devtools::load_all()
test_that("Testing rescale_extracted navios vs lx", {
    test_datadir <- file.path(testthat::test_path(), "testdata", "gated_lymphocytes")
    relevant_cols <- c(
        "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
        "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
    )
    # checking_devices_samples <- c("navios", "lx", "s001", "s051", "s132")
    checking_devices_samples <- c("navios", "lx", "s001", "s132")
    extracted_mfis <- sapply(checking_devices_samples, simplify = FALSE, function(device_x) {
        extract_mfi(
            fcs_dir = file.path(test_datadir, device_x),
            regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
            transform = function(x) {
                asinh(x / 5e3)
            },
            # gating_set_file = "gates.removeme",
            # gate_extract = "/Lymphocytes",
            multistaining = TRUE,
            regex_multistain = "(_15-MasterMix)\\.fcs",
            multistain_columns = relevant_cols
        )
    })


    extracted_quantile <- sapply(checking_devices_samples, simplify = FALSE, function(device_x) {
        extract_quantile(
            fcs_dir = file.path(test_datadir, device_x),
            regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
            transform = function(x) {
                asinh(x / 5e3)
            },
            # gating_set_file = "gates.removeme",
            # gate_extract = "/Lymphocytes",
            multistaining = TRUE,
            regex_multistain = "(_15-MasterMix)\\.fcs",
            multistain_columns = relevant_cols
        )
    })

    # # Rescale SINGLE NONE samples #####
    # loaded_singlenone <- list.files(
    #     path = test_datadir,
    #     # pattern = "(SingleNone)\\.fcs",
    #     # pattern = "(CD3-AA750)\\.fcs",
    #     # pattern = "(CD3-AF700)\\.fcs",
    #     pattern = "(CD3-FITC)\\.fcs",
    #     full.names = TRUE,
    #     recursive = TRUE
    # ) |>
    #     sapply(flowCore::read.FCS, simplify = FALSE)
    # names(loaded_singlenone) <- basename(dirname(names(loaded_singlenone)))
    # loaded_singlenone <- loaded_singlenone[checking_devices_samples]
    
    # rescaled_singlenone <- lapply(
    #     list(
    #         "raw" = scale_column_identity,
    #         "minmax" = scale_column_minmax,
    #         "relative" = scale_column_relative
    #     ),
    #     function(scale_column_fun) {
    #         sapply(
    #             names(loaded_singlenone),
    #             simplify = FALSE,
    #             function(device) {
    #                 rescale_extracted(
    #                     sample_to_rescale = flowCore::exprs(loaded_singlenone[[device]]),
    #                     extracted_mfi = extracted_mfis[[device]],
    #                     scale_column_fun = scale_column_fun
    #                 )
    #             }
    #         )
    #     }
    # )

    # # Quantile relativisation requires ONLY the cd3 positive signals
    # loaded_singlenone_cd3 <- lapply(
    #     loaded_singlenone,
    #     function(x) {
    #         clustering <- stats::kmeans(flowCore::exprs(x[, relevant_cols]), centers = 2)
    #         center_high <- rownames(clustering$centers)[which.max(apply(clustering$centers, 1, sum))]
    #         x[clustering$cluster == center_high, ]
    #     }
    # )
    # rescaled_singlenone_quantile <- lapply(
    #     list(
    #         "Q.minmax" = scale_column_minmaxQ,
    #         "Q.relative" = scale_column_relativeQ
    #     ),
    #     function(scale_column_fun) {
    #         sapply(
    #             names(loaded_singlenone_cd3),
    #             simplify = FALSE,
    #             function(device) {
    #                 tmp <- extracted_quantile[[device]]$single
    #                 names(tmp) <- NULL
    #                 rescaled_sample <- rescale_quantile(
    #                     sample_to_rescale = flowCore::exprs(loaded_singlenone_cd3[[device]]),
    #                     extracted_mfi = extracted_mfis[[device]],
    #                     extracted_qfun = unlist(tmp, recursive = FALSE),
    #                     scale_column_fun = scale_column_fun, 
    #                     # reference_marker = "AF700-A"
    #                     reference_marker = "FITC-A"
    #                 )
    #                 rescaled_sample
    #             }
    #         )
    #     }
    # )
    # rescaled_samples_singlenone <- c(
    #     rescaled_singlenone,
    #     rescaled_singlenone_quantile
    # )
    # pdf("removeme-singlenone.pdf", width = 25, height = 12)
    # print(
    #     plot_compare_alignment(rescaled_samples_singlenone, relevant_cols, ylim = NA, faceting = "vertical")
    # )
    # dev.off()


    # Rescale PANEL samples #####
    loaded_panel <- list.files(
        path = file.path(testthat::test_path(), "testdata", "gated_cd3"),
        pattern = "(panel)\\.fcs",
        full.names = TRUE,
        recursive = TRUE
    ) |>
        sapply(flowCore::read.FCS, simplify = FALSE)
    names(loaded_panel) <- basename(dirname(names(loaded_panel)))
    loaded_panel <- loaded_panel[checking_devices_samples]

    rescaled_panelsample <- lapply(
        list(
            "raw" = scale_column_identity,
            "minmax" = scale_column_minmax,
            "relative" = scale_column_relative
        ),
        function(scale_column_fun) {
            sapply(
                names(loaded_panel),
                simplify = FALSE,
                function(device) {
                    rescale_extracted(
                        sample_to_rescale = flowCore::exprs(loaded_panel[[device]]),
                        extracted_mfi = extracted_mfis[[device]],
                        scale_column_fun = scale_column_fun
                    )
                }
            )
        }
    )


    rescaled_panelsample_quantile <- lapply(
        list(
            "Q.minmax" = scale_column_minmaxQ,
            "Q.relative" = scale_column_relativeQ
        ),
        function(scale_column_fun) {
            sapply(
                names(loaded_panel),
                simplify = FALSE,
                function(device) {
                    pdf(paste0("removeme-", device, ".pdf"))
                    tmp <- extracted_quantile[[device]]$single
                    names(tmp) <- NULL
                    rescaled_sample <- rescale_quantile(
                        sample_to_rescale = flowCore::exprs(loaded_panel[[device]]),
                        extracted_mfi = extracted_mfis[[device]],
                        extracted_qfun = unlist(tmp, recursive = FALSE),
                        scale_column_fun = scale_column_fun
                    )
                    dev.off()
                    rescaled_sample
                }
            )
        }
    )

    rescaled_samples <- c(
        rescaled_panelsample,
        rescaled_panelsample_quantile
    )
    pdf("removeme.pdf", width = 25, height = 12)
    print(
        plot_compare_alignment(rescaled_samples, relevant_cols, ylim = NA, faceting = "vertical")
    )
    dev.off()
    browser()
})
