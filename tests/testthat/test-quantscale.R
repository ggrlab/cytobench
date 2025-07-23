# test_that("Testing rescale_extracted", {
#     extracted_mfis <- extract_singlestain_mfi(
#         "data-raw/s001",
#         gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
#     )
#     loaded_fcs <- flowCore::read.FCS("data-raw/s001/CT_p001_s001_comp-manual_ungated_none_100k_navios_14-01..11-SingleNone.fcs")
#     rescaled_sample <- rescale_extracted(
#         sample_to_rescale = flowCore::exprs(loaded_fcs),
#         extracted_mfi = extracted_mfis
#     )
#     testthat::expect(TRUE)
# })

devtools::load_all()
test_that("Testing rescale_extracted navios vs lx", {
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


#     extracted_mfis <- sapply(c("navios", "lx"), simplify = FALSE, function(device_x) {
#         extract_mfi(
#             fcs_dir = file.path(test_datadir, device_x),
#             regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
#             transform = function(x) {
#                 asinh(x / 5e3)
#             },
#             # gating_set_file = "gates.removeme",
#             # gate_extract = "/Lymphocytes",
#             multistaining = TRUE,
#             regex_multistain = "(_15-MasterMix)\\.fcs"
#         )
#     })


#     extracted_mfis_quantized <- sapply(c("navios", "lx"), simplify = FALSE, function(device_x) {
#         extract_mfi_quantized(
#             fcs_dir = file.path(test_datadir, device_x),
#             regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
#             transform = function(x) {
#                 asinh(x / 5e3)
#             },
#             # gating_set_file = "gates.removeme",
#             # gate_extract = "/Lymphocytes",
#             multistaining = TRUE,
#             regex_multistain = "(_15-MasterMix)\\.fcs"
#         )
#     })
#     stop()

#     # The reatios are just for information
#     ratio_mfis <- extracted_mfis[["lx"]][, c("negative", "positive", "negative.multi", "positive.multi")] /
#         extracted_mfis[["navios"]][, c("negative", "positive", "negative.multi", "positive.multi")]


#     loaded_singlenone <- list.files(
#         path = test_datadir,
#         pattern = "(SingleNone)\\.fcs",
#         full.names = TRUE,
#         recursive = TRUE
#     ) |>
#         sapply(flowCore::read.FCS, simplify = FALSE)
#     names(loaded_singlenone) <- basename(dirname(names(loaded_singlenone)))

#     rescaled_singlenone <- lapply(
#         list(
#             "raw" = scale_column_identity,
#             "minmax" = scale_column_minmax,
#             "relative" = scale_column_relative
#         ),
#         function(scale_column_fun) {
#             sapply(
#                 names(loaded_singlenone),
#                 simplify = FALSE,
#                 function(device) {
#                     rescale_extracted(
#                         sample_to_rescale = flowCore::exprs(loaded_singlenone[[device]]),
#                         extracted_mfi = extracted_mfis[[device]],
#                         scale_column_fun = scale_column_fun
#                     )
#                 }
#             )
#         }
#     )
#     rs_sn_long <- lapply(rescaled_singlenone, data.table::rbindlist, idcol = "device", fill = TRUE)
#     rescaled_variants <- data.table::rbindlist(rs_sn_long, idcol = "scale_column_fun", fill = TRUE)
#     rescaled_variants_long <- data.table::melt(
#         rescaled_variants,
#         id.vars = c("device", "scale_column_fun"),
#         variable.name = "feature",
#         value.name = "value"
#     )
#     densfuns_ranges <- rescaled_variants_long[,
#         {
#             if (scale_column_fun == "raw") {
#                 value <- asinh(value / 1e3)
#             } else {
#                 value <- asinh(value * 100)
#             }
#             resfun <- function(x) {
#                 NA
#             }
#             try(resfun <- approxfun(density(value, na.rm = TRUE)), silent = FALSE)
#             .(
#                 dens_approxfun = list(list(resfun)),
#                 min = min(value),
#                 max = max(value)
#             )
#         },
#         by = c("device", "scale_column_fun", "feature")
#     ]
#     n_resolution <- 1e3
#     approximated_dens_all <- data.table()
#     for (feature_x in unique(densfuns_ranges$feature)) {
#         current <- densfuns_ranges[feature == feature_x]
#         range_total <- range(current$min, current$max, na.rm = TRUE)

#         approximated_dens <- data.table()
#         for (i in seq_along(current$dens_approxfun)) {
#             afun <- current[["dens_approxfun"]][[i]][[1]]
#             range_current <- c(range_total[1], range_total[2])
#             x <- seq(range_current[[1]], range_current[[2]], length.out = n_resolution)
#             approximated_dens <- rbind(
#                 approximated_dens, cbind(
#                     current[i],
#                     data.table(
#                         x = x,
#                         y = afun(x)
#                     )
#                 )
#             )
#         }
#         approximated_dens_all <- rbind(
#             approximated_dens_all,
#             approximated_dens
#         )
#     }
#     approximated_dens_all[, scale_column_fun := factor(scale_column_fun, levels = unique(c("raw", "minmax", "relative", names(rescaled_singlenone))))]

#     #  [1] "FSC-A"   "FSC-W"   "SSC-H"   "SSC-A"   "FITC-A"  "PE-A"    "ECD-A"
#     #  [8] "PC5.5-A" "PC7-A"   "APC-A"   "AF700-A" "AA750-A" "PB-A"    "KrO-A"
#     # [15] "TIME"    "FSC-H"
#     relevant_features <- c(
#         "FITC-A", "PE-A", "ECD-A", "PC5.5-A",
#         "PC7-A", "APC-A", "AF700-A", "AA750-A",
#         "PB-A", "KrO-A"
#     )
#     pdf("removeme.pdf", width = 15, height = 25)
#     print(
#         ggplot2::ggplot(
#             approximated_dens_all[feature %in% relevant_features], ggplot2::aes(x = x, y = y, col = device)
#         ) +
#             ggplot2::geom_line() +
#             ggplot2::ylim(c(0, .4)) +
#             # ggplot2::geom_area(alpha = .5) +
#             # ggplot2::facet_wrap(scale_column_fun~feature, scales = "free") +
#             # ggplot2::facet_grid(feature ~ scale_column_fun, scales = "free_y") +
#             ggh4x::facet_grid2(
#                 feature ~ scale_column_fun,
#                 scales = "free",
#                 independent = "all"
#             ) +
#             ggpubr::theme_pubclean()
#     )
#     dev.off()
# })
