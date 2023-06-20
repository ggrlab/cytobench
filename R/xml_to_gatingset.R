xml_to_gatingset <- function(read_xml, analysis_fcs_file, verbose = FALSE, verbose_pdf_path = "removeme_extracted_gating.pdf") {
    # Extract sample-by-sample:
    # //  finds anywhere in the document (ignoring the current node)
    # .// finds anywhere beneath the current node
    all_entries <- xml2::xml_find_all(read_xml, xpath = ".//Entries")
    datasets <- xml2::xml_find_all(all_entries, xpath = ".//DataSet")
    dataset_names <- xml2::xml_attr(datasets, "N")
    gates <- xml2::xml_find_all(all_entries, xpath = ".//Gates")
    gate_constraints <- xml2::xml_find_all(all_entries, xpath = ".//GateConstraints")

    gatelist_per_sample <- lapply(gates, xmlgates_to_gatelist, gate_constraints = gate_constraints)
    names(gatelist_per_sample) <- dataset_names

    # take one random .fcs, necessary to create a GatingSet
    analysis_fcs <- flowCore::read.FCS(
        analysis_fcs_file,
        truncate_max_range = FALSE,
        # Read only a single cell
        which.lines = 1
    )
    zero_cell_flowframe <- analysis_fcs[0, ]

    gatingset_per_sample <- lapply(
        gatelist_per_sample,
        gatelist_to_gatingset,
        flowframe = zero_cell_flowframe
    )

    if (verbose) {
        if (length(gatingset_per_sample) > 1) {
            warning("Verboseness is only showing the first of more than 1 samples!")
        }
        print(flowWorkspace::gs_get_pop_paths(gatingset_per_sample[[1]], path = "full"))
        library(flowWorkspace)
        pdf(verbose_pdf_path)
        flowWorkspace::plot(gatingset_per_sample[[1]])
        dev.off()
        cat("Wrote", verbose_pdf_path, "\n")
    }
    return(gatingset_per_sample)
}

gatelist_to_gatingset <- function(gatelist, flowframe) {
    zero_gs <- flowWorkspace::GatingSet(flowCore::flowSet(flowframe))
    i <- 1
    while (length(gatelist) >= 1) {
        try({
            flowWorkspace::gs_pop_add(
                zero_gs,
                gatelist[[i]][["gate"]],
                parent = ifelse(is.null(gatelist[[i]][["parent"]]), "root", gatelist[[i]][["parent"]])
            )
            gatelist[[i]] <- NULL
            i <- i - 1
        })
        i <- i + 1
        if (i > length(gatelist)) {
            i <- 1
        }
    }
    return(zero_gs)
}

xmlgates_to_gatelist <- function(xml_gates_single, gate_constraints) {
    # Interpret the gates
    # gate_contents <- xml2::xml_contents(xml_gates_single)
    # single_gate <- xml2::xml_contents(gate_contents[[1]])
    gates_list <- xml2::as_list(xml_gates_single)
    scale_interpretation <- c("I" = "linear", "O" = "log", "C" = "logicle")
    interpreted_gates <- lapply(gates_list, function(single_gate) {
        markers <- lapply(single_gate[["A"]], function(x) {
            c("name" = attributes(x)[["M"]], "scale" = scale_interpretation[attributes(x)[["S"]]])
        })

        coordinates <- lapply(single_gate[["P"]], function(x) {
            x_y <- as.numeric(strsplit(attributes(x)[["O"]], ",")[[1]])
            c("x" = x_y[1], "y" = x_y[2])
        })

        radius <- list()
        if ("Coefficients" %in% names(single_gate)) {
            radius <- lapply(single_gate[["Coefficients"]], function(x) {
                x_y <- as.numeric(x)
                return(x_y)
            })
            names(radius) <- c("x", "y")
        }
        return(
            list("markers" = markers, "coordinates" = coordinates, "radius" = radius, "attr" = attributes(single_gate))
        )
    })

    gate_names <- unlist(lapply(interpreted_gates, function(x) x[["attr"]][["N"]]))
    if (anyDuplicated(gate_names)) {
        warning("Duplicated gate names, setting new names")
        gate_names <- paste0(gate_names, "_", seq_len(length(gate_names)))
    }
    # L: one-way linear
    # P: Polygon
    # R: Rectangle
    # E: Ellipse
    gate_constraints_content <- xml2::xml_contents(gate_constraints)
    gate_constraints_l <- xml2::as_list(gate_constraints_content)
    names(gate_constraints_l) <- xml2::xml_name(gate_constraints_content)
    histogram_dividers <- gate_constraints_l[names(gate_constraints_l) == "HistogramDivider"]
    histogram_divider_gates_right <- lapply(histogram_dividers, function(x) attr(x, "RG"))
    histogram_divider_gates_left <- lapply(histogram_dividers, function(x) attr(x, "LG"))

    staggered <- gate_constraints_l[names(gate_constraints_l) == "Staggered"]
    staggered_grouped <- lapply(c("topleft" = "ULG", "topright" = "URG", "bottomleft" = "LLG", "bottomright" = "LRG"), function(x) {
        lapply(staggered, function(y) attr(y, x))
    })
    staggered_grouped_unlist <- unlist(staggered_grouped)
    # "Staggered" seem to be the Quadrant-gates

    # If $attr contains "G" as attr(ibute), it is a child
    # $attr$N is the label of the gate
    gatelist <- list()
    for (gate_i in seq_len(length(interpreted_gates))) {
        gate_type <- names(interpreted_gates[gate_i])
        single_gate <- interpreted_gates[[gate_i]]
        single_gate_list <- lapply(single_gate[["markers"]], function(x) x[["name"]])
        names(single_gate_list) <- unlist(single_gate_list)

        coords <- do.call(rbind, single_gate[["coordinates"]])
        if (gate_type == "L") {
            part_coord <- coords[, "x", drop = TRUE]
            if (gate_names[gate_i] %in% histogram_divider_gates_left) {
                part_coord[which.min(part_coord)] <- -Inf
                # print("LEFT gate")
            }
            if (gate_names[gate_i] %in% histogram_divider_gates_right) {
                part_coord[which.max(part_coord)] <- Inf
                # print("RIGHT gate")
            }
            single_gate_list[[1]] <- part_coord
            fc_gate <- flowCore::rectangleGate(filterId = gate_names[gate_i], single_gate_list)
        } else if (gate_type == "P") {
            # if (gate_names[gate_i] %in% staggered_grouped_unlist) {
            #     # I started this here but it is not finished
            #     current_name <- names(staggered_grouped_unlist)[staggered_grouped_unlist == gate_names[gate_i]]
            #     if (startsWith("topright", current_name)) {
            #         coords[coords[, "x"] == max(coords[, "x"]), "x"] <- Inf
            #     } else if (startsWith("bottomright", current_name)) {
            #         coords[coords[, "x"] == max(coords[, "x"]), "x"] <- Inf
            #     } else if (startsWith("topleft", current_name)) {
            #         coords[coords[, "x"] == min(coords[, "x"]), "x"] <- -Inf
            #     } else if (startsWith("bottomleft", current_name)) {
            #         coords[coords[, "x"] == min(coords[, "x"]), "x"] <- -Inf
            #         coords[coords[, "y"] == min(coords[, "y"]), "y"] <- -Inf
            #     }
            #     # print("LEFT gate")
            # }

            colnames(coords) <- names(single_gate_list)
            fc_gate <- flowCore::polygonGate(filterId = gate_names[gate_i], coords)
            # Warning: The 'boundaries' argument is deprecated, please use '.gate' instead.
            # # from ?flowCore::polygonGate
            # ## Defining the gate
            # sqrcut <- matrix(c(300,300,600,600,50,300,300,50),ncol=2,nrow=4)
            # colnames(sqrcut) <- c("FSC-H","SSC-H")
            # pg <- polygonGate(filterId="nonDebris", boundaries= sqrcut)
            # pg
            if (all(coords == 0)) {
                warning(paste0(
                    "Gate ",
                    gate_names[gate_i],
                    ": Hinged gate is giving wrong sizes. I do not know what Kaluza is exactly ",
                    "doing. To press forward now: Do not use hinged gates, no cell will be gated but the gate created "
                ))
            }
        } else if (gate_type == "R") {
            colnames(coords) <- names(single_gate_list)
            fc_gate <- flowCore::rectangleGate(filterId = gate_names[gate_i], coords)
        } else if (gate_type == "E") {
            warning(paste0(
                "Gate ",
                gate_names[gate_i],
                ": Ellipsoid gate is still giving wrong sizes. I do not know what Kaluza is exactly ",
                "doing. To press forward now: Rescaling any axis of the ellipse renders it as ",
                "many-point-polygon within Kaluza which I can read easily. Do this in the end for ",
                "all your ellipses."
            ))
            # https://support.bioconductor.org/p/35360/
            cov.matrix <- function(a, b, angle) {
                theta <- angle * (pi / 180)
                c1 <- ((cos(theta)^2) / a^2) + ((sin(theta)^2) / b^2)
                c2 <- sin(theta) * cos(theta) * ((1 / a^2) - (1 / b^2))
                c3 <- ((sin(theta)^2) / a^2) + ((cos(theta)^2) / b^2)

                m1 <- matrix(c(c1, c2, c2, c3), byrow = TRUE, ncol = 2)
                # print(m1)
                m2 <- solve(m1)

                m2
            }
            coords <- single_gate[["coordinates"]][["P"]]
            radius <- unlist(single_gate[["radius"]])
            # The following number is completely arbitrary, I have absolutely no idea where
            # it comes from... was empirically tested by
            # 1. drawing a kaluza gate and counting the cells inside
            # 2. Taking the circle-gate then apply it by myself
            # 3. move it around until it fits.
            radius <- radius / 144.65
            # radius <- radius / 100
            # angle <- as.numeric(single_gate$attr[["L"]]) + 90
            angle <- as.numeric(single_gate$attr$L)
            if (length(angle) == 0) {
                angle <- 0
            }
            # angle_radian <- (angle * pi) / 180

            # sigma_x <- radius["x"] / sqrt(2)
            # sigma_y <- radius["y"] / sqrt(2)
            # rho <- cos(angle_radian) * sin(angle_radian)

            # covariance_matrix <- matrix(
            #     c(
            #         sigma_x^2, rho * sigma_x * sigma_y,
            #         rho * sigma_x * sigma_y, sigma_y^2
            #     ),
            #     byrow = TRUE, ncol = 2
            # )
            covariance_matrix <- cov.matrix(
                a = radius["x"],
                b = radius["y"],
                angle = angle
            )
            colnames(covariance_matrix) <- names(single_gate_list)
            rownames(covariance_matrix) <- names(single_gate_list)
            fc_gate <- flowCore::ellipsoidGate(filterId = gate_names[gate_i], covariance_matrix, mean = c(coords))
        } else if (gate_type == "B") {
            # Boolean gates
            stop("Boolean gates are not implemented yet.")
        } else {
            print(interpreted_gates[gate_i])
            print(single_gate)
            stop(paste0("Unrecognized gate type '", gate_type, "'"))
            # stop("Unrecognized gate type")
        }
        gatelist[[gate_i]] <- list("gate" = fc_gate, "parent" = single_gate$attr[["G"]])
    }
    return(gatelist)
}
