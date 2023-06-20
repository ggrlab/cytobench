#' @export
kaluza_check_gated <- function(obj_kaluza_read_analysis,
                               exported_gates_obj_or_path,
                               gating_pdf = "exported_gating.pdf",
                               add_gating_tables = TRUE) {
    pop <- Number <- count <- NULL # Avoid lint warnings
    if (length(obj_kaluza_read_analysis[["applied_gatings"]]) == 0) {
        stop(paste0(
            "Please set apply_fcs = TRUE in kaluza_read_analysis.",
            " You want to compare the applied gating counts to the true ones. "
        ))
    }
    if (length(exported_gates_obj_or_path) == 1 && file.exists(exported_gates_obj_or_path)) {
        kaluza_true_counts <- tibble::as_tibble(data.table::fread(
            exported_gates_obj_or_path,
            sep = ";", dec = ",",
            skip = 3
        ))
    } else {
        kaluza_true_counts <- exported_gates_obj_or_path
    }
    if (anyDuplicated(kaluza_true_counts$Gate) != 0) {
        stop("Duplicated gate names in example_fcs Ungated Statistics.csv - cannot merge according to them")
    }


    extracted_counts <- tibble::as_tibble(
        flowWorkspace::gs_pop_get_stats(obj_kaluza_read_analysis[["applied_gatings"]])
    )

    extracted_counts <- extracted_counts %>% dplyr::mutate(
        Gate = stringr::str_replace(pop, ".*/", "")
    )
    extracted_counts$Gate[extracted_counts$Gate == "root"] <- "All"

    joined_true_mine <- dplyr::full_join(kaluza_true_counts, extracted_counts, by = "Gate")
    different_gates <- joined_true_mine %>%
        dplyr::filter(
            # Number comes from kaluza_true_counts (and is therefore "true")
            # Count comes from extracted_counts (and is therefore "mine")
            Number != count
        ) %>%
        dplyr::mutate(ratio_true_mine = Number / count)

    if (require("ggcyto")) {
        gh <- obj_kaluza_read_analysis[["applied_gatings"]][[1]][[1]]
        nodes <- flowWorkspace::gs_get_pop_paths(gh)[-1] # exclude root node

        pdf(gating_pdf)
        for (node_x in nodes) {
            cat("\n\n")
            p1 <- ggcyto::autoplot(gh, node_x) +
                # clear all the geom_stats() layer previously added
                ggcyto::stats_null() +
                ggcyto::scale_x_flowCore_fasinh() +
                ggcyto::scale_y_flowCore_fasinh() +
                ggcyto::geom_stats(type = c("gate_name", "count"))
            print(p1)
            if (require("ggpp") && add_gating_tables) {
                one_gate <- flowWorkspace::gs_pop_get_gate(gh, node_x)[[1]]
                if ("rectangleGate" %in% class(one_gate)) {
                    boundaries <- rbind(one_gate@min, one_gate@max)
                    boundaries <- tibble::as_tibble(boundaries)
                } else if ("polygonGate" %in% class(one_gate)) {
                    boundaries <- one_gate@boundaries
                    boundaries <- tibble::as_tibble(boundaries)
                } else {
                    warning("Gate cannot be tabled: ", str(one_gate))
                    next
                }
                df <- tibble::tibble(x = 0, y = 0, tb = list(boundaries))

                p2 <- ggplot2::ggplot() +
                    ggpp::geom_table(
                        data = df,
                        ggplot2::aes(x = x, y = y, label = tb)
                    ) +
                    ggplot2::theme_void()
                print(p2)
            }
        }
        dev.off()
        cat("Wrote ", gating_pdf, "\n")
    }
    if (nrow(different_gates) > 0) {
        print(different_gates, n = 10000)
        stop("The previous gates did not match. Potentially the extracted gating was not correct.")
    }
    return(TRUE)
}
