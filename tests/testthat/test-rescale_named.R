test_that("Testing rescale_extracted", {
    set.seed(42)
    fs_ss <- simulate_cd3()
    tmpdir <- withr::local_tempdir()
    flowCore::write.flowSet(fs_ss, tmpdir)

    extracted_mfis_singlestain <- extract_mfi(
        tmpdir,
        regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$"
    )

    rescaled_sample <- rescale_extracted(
        sample_to_rescale = flowCore::exprs(fs_ss[["sample0_12-panel"]]),
        extracted_mfi = extracted_mfis_singlestain
    )
    testthat::expect(TRUE)
})

test_that("Testing rescale_extracted", {
    browser()
    set.seed(42)
    fs_ss <- simulate_cd3()
    set.seed(43)
    fs_ss_2 <- simulate_cd3()
    tmpdir <- withr::local_tempdir()
    flowCore::write.flowSet(fs_ss, tmpdir)

    extracted_mfis <- extract_mfi(
        tmpdir,
        regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$"
    )
    # loaded_fcs <- flowCore::read.FCS("data-raw/s001/CT_p001_s001_comp-manual_ungated_none_100k_navios_14-01..11-SingleNone.fcs")
    loaded_cy1 <- fs_ss[["sample0_12-panel"]]
    loaded_cy2 <- fs_ss_2[["sample0_12-panel"]]


    tmpplot <- function(df_x, asdf) {
        ggplot2::ggplot(
            df_x, ggplot2::aes(x = !!ggplot2::sym(asdf))
        ) +
            ggplot2::geom_density(fill = "blue", alpha = .3) +
            ggplot2::ylim(c(0, .1)) +
            ggpubr::theme_pubclean()
    }

    pc7_neg <- extracted_mfis[extracted_mfis$feature == "PC7-A", "negative"][[1]]
    pc7_pos <- extracted_mfis[extracted_mfis$feature == "PC7-A", "positive"][[1]]

    plot(((1:1000) / (10)))
    x <- (-1000:1000) / 200
    plot(x, (exp(x) / (1 + exp(x))))
    absolute_sigmoid <- function(x) {
        ((1 / (1 + exp(-x))) - .5) * 2
    }
    plot(x, absolute_sigmoid(x))


    refactor <- 1e3
    df <- tibble::tibble(
        x = c(
            flowCore::exprs(loaded_cy1)[, "PC7-A"],
            flowCore::exprs(loaded_cy2)[, "PC7-A"]
        ),
        device = c(rep("lx", nrow(loaded_cy1)), rep("navios", nrow(loaded_cy2)))
    ) |>
        dplyr::mutate(
            x_asinh = asinh(x / 500),
            x_minmax = (x - min(x)) / (max(x) - min(x)),
            x_minmax_asinh = asinh(x_minmax * refactor),
            x_aligned = (x - pc7_neg) / (pc7_pos - pc7_neg),
            x_aligned_asinh = asinh(x_aligned * refactor),
            x_aligned_v2 = ((x - pc7_neg) / (pc7_pos - pc7_neg)) * abs(absolute_sigmoid(x)) * refactor,
            x_aligned_v2_asinh = asinh(x_aligned_v2),
            x_aligned_v3_asinh = asinh(absolute_sigmoid(x_aligned) * refactor)
        )
    tmpdir <- withr::local_tempdir()
    pdf(file.path(tmpdir, "removeme.pdf"))
    print(tmpplot(df, "x_asinh"))
    print(tmpplot(df, "x_minmax_asinh"))
    print(tmpplot(df, "x_aligned_asinh"))
    print(tmpplot(df, "x_aligned_v2_asinh"))
    print(tmpplot(df, "x_aligned_v3_asinh"))
    dev.off()

    ##### PB:
    pb_neg <- extracted_mfis[extracted_mfis$feature == "PB-A", "negative"][[1]]
    pb_pos <- extracted_mfis[extracted_mfis$feature == "PB-A", "positive"][[1]]

    plot(((1:1000) / (10)))
    x <- (-1000:1000) / 200
    plot(x, (exp(x) / (1 + exp(x))))
    absolute_sigmoid <- function(x) {
        ((1 / (1 + exp(-x))) - .5) * 2
    }
    plot(x, absolute_sigmoid(x))


    refactor <- 1e2
    df <- tibble::tibble(x = flowCore::exprs(loaded_cy1)[, "PB-A"]) |>
        dplyr::mutate(
            x_asinh = asinh(x / 500),
            x_minmax = (x - min(x)) / (max(x) - min(x)),
            x_minmax_asinh = asinh(x_minmax * refactor),
            x_aligned = (x - pb_neg) / (pb_pos - pb_neg),
            x_aligned_asinh = asinh(x_aligned * refactor),
            x_aligned_v2 = ((x - pb_neg) / (pb_pos - pb_neg)) * abs(absolute_sigmoid(x)) * refactor,
            x_aligned_v2_asinh = asinh(x_aligned_v2),
            x_aligned_v3_asinh = asinh(absolute_sigmoid(x_aligned) * refactor)
        )
    tmpdir <- withr::local_tempdir()
    pdf(file.path(tmpdir, "removeme.pdf"))
    print(tmpplot(df, "x_asinh"))
    print(tmpplot(df, "x_minmax_asinh"))
    print(tmpplot(df, "x_aligned_asinh"))
    print(tmpplot(df, "x_aligned_v2_asinh"))
    print(tmpplot(df, "x_aligned_v3_asinh"))
    dev.off()
    rescaled_sample <- rescale_extracted(
        sample_to_rescale = flowCore::exprs(loaded_cy1),
        extracted_mfi = extracted_mfis
    )
    testthat::expect(TRUE)
})

#### I have removed testdata/gated_lymphocytes as the data was too large for CRAN

# test_that("Testing rescale_extracted navios vs lx", {
#     test_datadir <- file.path(testthat::test_path(), "testdata", "gated_lymphocytes")
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
#             "minmax" = scale_column_minmax,
#             "minmax_asinh" = scale_column_minmax_asinh,
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
#             value <- asinh(value * 1e2)
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


#     #  [1] "FSC-A"   "FSC-W"   "SSC-H"   "SSC-A"   "FITC-A"  "PE-A"    "ECD-A"
#     #  [8] "PC5.5-A" "PC7-A"   "APC-A"   "AF700-A" "AA750-A" "PB-A"    "KrO-A"
#     # [15] "TIME"    "FSC-H"
#     relevant_features <- c(
#         "FITC-A", "PE-A", "ECD-A", "PC5.5-A",
#         "PC7-A", "APC-A", "AF700-A", "AA750-A",
#         "PB-A", "KrO-A"
#     )
#     pdf("removeme.pdf", width = 25, height = 25)
#     ggplot2::ggplot(
#         approximated_dens_all[feature %in% relevant_features], ggplot2::aes(x = x, y = y, col = device)
#     ) +
#         ggplot2::geom_line() +
#         # ggplot2::facet_wrap(scale_column_fun~feature, scales = "free") +
#         # ggplot2::facet_grid(feature ~ scale_column_fun, scales = "free_y") +
#         ggh4x::facet_grid2(
#             feature ~ scale_column_fun,
#             scales = "free",
#             independent = "all"
#         ) +
#         ggpubr::theme_pubclean()
#     dev.off()
# })
