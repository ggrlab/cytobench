#' @export 
kaluza_read_analysis <- function(path, apply_fcs = FALSE, verbose = TRUE, which_xml = 0) {
    current_tempdir <- tempdir_time()
    if (length(list.files(current_tempdir)) > 0) {
        stop("Tempdir is not empty, try again, my tempdir is including the time.")
    }
    ## 1. unzip into a temporary directory
    if (verbose) {
        cat("------ Unzipping -----\n")
    }
    unzip(path, exdir = current_tempdir)
    # results in
    #   - .fcs files for every dataset in Kaluza
    #   - One .xml
    #   - One .app

    ## 2. Extract the xml file
    if (verbose) {
        cat("------ Extracting .XML -----\n")
    }
    analysis_xml <- list.files(current_tempdir, pattern = "*.xml$", full.names = TRUE)
    if (length(analysis_xml) > 1 && which_xml == 0) {
        stop(paste0(
            "More than 1 xml file (",
            length(analysis_xml),
            ") unsupported, give the number which_xml should be used. ",
            "VERY cautious, I don't know how that could have happened."
        ))
    } else {
        which_xml <- 1
    }
    read_xml <- xml2::read_xml(analysis_xml[which_xml])

    if (verbose) {
        xml2::write_xml(read_xml, "verbose_xml.xml")
        cat("Wrote ./verbose_xml.xml\n")
    }

    all_fcs <- list.files(current_tempdir, pattern = "*.fcs", full.names = TRUE)
    if (verbose) {
        cat("------ xml to gatingset -----\n")
    }
    gs <- xml_to_gatingset(
        read_xml = read_xml,
        analysis_fcs_file = all_fcs[1], # take any fcs
        verbose = verbose
    )

    applied_gating_list <- list()
    if (apply_fcs) {
        if (verbose) {
            cat("------ Reading + Applying gatingsets to all fcs -----\n")
        }
        for (sample_i in seq_len(length(gs))) {
            read_fcs <- flowWorkspace::load_cytoset_from_fcs(all_fcs[[sample_i]])
            applied_gating_list[[sample_i]] <- suppressMessages(
                flowWorkspace::gh_apply_to_cs(
                    x = gs[[sample_i]][[1]],
                    cs = read_fcs,
                )
            )
        }
        names(applied_gating_list) <- names(gs)
    }
    return(list(
        "gatingsets" = gs,
        "applied_gatings" = applied_gating_list
    ))
}
