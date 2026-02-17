
require_namespace <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
    }
}

validate_pipeline_spec <- function(spec) {
    required_top <- c("name", "input", "output", "gating", "plotting")
    missing_top <- setdiff(required_top, names(spec))
    if (length(missing_top) > 0) {
        stop(sprintf("Missing top-level spec entries: %s", paste(missing_top, collapse = ", ")), call. = FALSE)
    }

    required_input <- c("analysis_paths", "read_protocol", "read_ds2")
    missing_input <- setdiff(required_input, names(spec$input))
    if (length(missing_input) > 0) {
        stop(sprintf("Missing input spec entries: %s", paste(missing_input, collapse = ", ")), call. = FALSE)
    }

    if (!is.function(spec$input$analysis_paths)) {
        stop("spec$input$analysis_paths must be a function returning analysis file paths.", call. = FALSE)
    }
    if (!is.function(spec$input$read_protocol)) {
        stop("spec$input$read_protocol must be a function(path) -> protocol.", call. = FALSE)
    }
    if (!is.function(spec$input$read_ds2)) {
        stop("spec$input$read_ds2 must be a function(path) -> named list of flowFrames.", call. = FALSE)
    }

    if (!is.function(spec$gating$relevant_gate)) {
        stop("spec$gating$relevant_gate must be a function(analysis_name, analysis_path) -> gate string.", call. = FALSE)
    }
    if (!is.function(spec$gating$gate_cells)) {
        stop("spec$gating$gate_cells must be a function(flowset, protocol, relevant_gate, spec, analysis_name) -> list.", call. = FALSE)
    }

    if (!is.function(spec$output$analysis_name)) {
        stop("spec$output$analysis_name must be a function(analysis_path) -> analysis_name.", call. = FALSE)
    }
    if (!is.function(spec$plotting$plot_fn)) {
        stop("spec$plotting$plot_fn must be a function(gated, relevant_gate, spec, analysis_path) -> list of ggplots.", call. = FALSE)
    }

    invisible(TRUE)
}

rename_by_ds1_position <- function(ff, ds1_reference, normalize_pattern = " ((LIN)|(LOG))", append_time_marker = TRUE) {
    require_namespace("flowCore")

    ref_markers <- flowCore::markernames(ds1_reference)
    if (append_time_marker) {
        ref_markers <- c(ref_markers, "TIME")
    }

    flowCore::markernames(ff) <- structure(unname(ref_markers), names = flowCore::colnames(ff))
    flowCore::colnames(ff) <- sub(normalize_pattern, "", flowCore::colnames(ds1_reference))
    ff
}

build_name_map <- function(ds1_reference, ds2_original, ds2_renamed) {
    require_namespace("flowCore")

    cbind(
        data.frame(
            colnames_ds1 = flowCore::colnames(ds1_reference),
            colnames_ds2 = flowCore::colnames(ds2_renamed),
            colnames_ds2_original = flowCore::colnames(ds2_original)
        ),
        rbind(
            data.frame(
                markernames_ds1 = flowCore::markernames(ds1_reference),
                markernames_ds2 = flowCore::markernames(ds2_renamed)
            ),
            NA
        )
    )
}

run_pipeline <- function(spec) {
    require_namespace("flowCore")
    require_namespace("flowWorkspace")
    require_namespace("data.table")

    validate_pipeline_spec(spec)

    analyses <- spec$input$analysis_paths()
    protocols <- stats::setNames(lapply(analyses, spec$input$read_protocol), analyses)

    for (analysis_path in analyses) {
        analysis_name <- spec$output$analysis_name(analysis_path)
        outdir <- file.path(spec$output$base_dir, analysis_name)
        dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)
        dir.create(file.path(outdir, "fcs"), showWarnings = FALSE, recursive = TRUE)

        cat(sprintf("Processing %s (%s)\n", analysis_name, basename(analysis_path)))

        ds1 <- NULL
        if (!is.null(spec$input$read_ds1) && is.function(spec$input$read_ds1)) {
            ds1 <- spec$input$read_ds1(analysis_path)
        }
        ds2 <- spec$input$read_ds2(analysis_path)

        ds2_compensated <- ds2
        if (!is.null(spec$compensation) && isTRUE(spec$compensation$enabled)) {
            ds2_compensated <- lapply(names(ds2), function(sample_x) {
                spec$compensation$apply(
                    flowframe = ds2[[sample_x]],
                    sample_name = sample_x,
                    protocol = protocols[[analysis_path]],
                    analysis_name = analysis_name,
                    analysis_path = analysis_path
                )
            })
            names(ds2_compensated) <- names(ds2)
        }

        ds2_renamed <- ds2_compensated
        name_maps <- NULL
        if (!is.null(spec$rename) && isTRUE(spec$rename$enabled)) {
            if (is.null(ds1)) {
                stop("Renaming is enabled but spec$input$read_ds1 is missing.", call. = FALSE)
            }
            ds1_reference <- ds1[[1]]
            ds2_renamed <- lapply(ds2_compensated, function(ff) {
                spec$rename$apply(ff = ff, ds1_reference = ds1_reference, spec = spec, analysis_name = analysis_name)
            })
            name_maps <- build_name_map(
                ds1_reference = ds1_reference,
                ds2_original = ds2[[1]],
                ds2_renamed = ds2_renamed[[1]]
            )
        }

        if (!is.null(spec$output$save_compensated) && isTRUE(spec$output$save_compensated)) {
            for (sample_name in names(ds2_renamed)) {
                outfile <- file.path(spec$output$compensated_dir, analysis_name, paste0(sample_name, ".fcs"))
                dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
                flowCore::write.FCS(ds2_renamed[[sample_name]], outfile)
            }
        }

        relevant_gate <- spec$gating$relevant_gate(analysis_name = analysis_name, analysis_path = analysis_path)
        gated <- spec$gating$gate_cells(
            flowset = flowCore::flowSet(ds2_renamed),
            protocol = protocols[[analysis_path]],
            relevant_gate = relevant_gate,
            spec = spec,
            analysis_name = analysis_name
        )

        gated_ff_list <- lapply(gated$flowset_gated, function(x) {
            tmp <- flowWorkspace::gs_pop_get_data(x, relevant_gate)
            tmp <- flowWorkspace::realize_view(tmp)
            flowWorkspace::cytoframe_to_flowFrame(tmp[[1]])
        })

        for (sample_name in names(gated_ff_list)) {
            flowCore::write.FCS(gated_ff_list[[sample_name]], file.path(outdir, "fcs", paste0(sample_name, ".fcs")))
        }

        if (!is.null(name_maps)) {
            data.table::fwrite(name_maps, file.path(outdir, "name_maps.csv"))
        }

        if (isTRUE(spec$plotting$enabled)) {
            plots <- spec$plotting$plot_fn(gated = gated, relevant_gate = relevant_gate, spec = spec, analysis_path = analysis_path)
            pdffile <- file.path(outdir, "plots", paste0(basename(analysis_path), "--gating_overview.pdf"))
            grDevices::pdf(pdffile, width = spec$plotting$pdf_width, height = spec$plotting$pdf_height)
            print(plots)
            grDevices::dev.off()
        }

        if (!is.null(spec$hooks$after_analysis) && is.function(spec$hooks$after_analysis)) {
            spec$hooks$after_analysis(
                analysis_path = analysis_path,
                analysis_name = analysis_name,
                gated = gated,
                gated_ff_list = gated_ff_list,
                ds2_renamed = ds2_renamed,
                outdir = outdir,
                spec = spec
            )
        }
    }

    invisible(TRUE)
}
