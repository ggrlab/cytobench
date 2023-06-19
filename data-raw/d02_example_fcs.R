## code to prepare `example_fcs` dataset goes here

example_csv <- data.table::fread(file.path("data-raw", "example_flowcyto.csv"), sep = ";", dec = ",")
# ff <- flowCore::flowFrame(as.matrix(example_csv))

# devtools::load_all()
export_fcs(
    list("example_fcs" = example_csv),
    safety_scaling = 1, safety_shift = 0,
    outdir = file.path("inst", "extdata"),
    new_colnames_to = "description",
    new_colnames = c(
        "FSC-A" = "no",
        "SSC-A" = "no",
        "FL1-A" = "marker",
        "FL4-A" = "marker",
        "FL5-A" = "marker",
        "FL7-A" = "marker",
        "FL8-A" = "marker",
        "FL9-A" = "marker",
        "FL10-A" = "marker",
        "FL11-A" = "marker",
        "FL12-A" = "marker",
        "FL13-A" = "marker",
        "FL14-A" = "marker",
        "FL16-A" = "marker",
        "FL21-A" = "marker",
        "FSC-Width" = "no",
        "TIME" = "no"
    )
)
