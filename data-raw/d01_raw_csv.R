## code to prepare `raw_csv` dataset goes here

paths_csvs <- list.files("data-raw/PRE_ARCSINH_CSV", full.names = TRUE, recursive = TRUE, pattern = "\\.csv$")
raw_csv <- lapply(paths_csvs, data.table::fread, sep = ";", dec = ",")
names(raw_csv) <- paths_csvs
set.seed(1209)
subset_csvs <- lapply(raw_csv, function(x) {
    x[sample(1:nrow(x), 500), ]
})
for (fp_x in names(raw_csv)) {
    fp_new <- gsub("data-raw/PRE_ARCSINH_CSV/", "inst/extdata/raw_csv/", fp_x)
    dir.create(dirname(fp_new), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(subset_csvs[[fp_x]], fp_new)
    cat("Wrote ", fp_new, "\n")
}
