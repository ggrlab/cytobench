#' @title Read All Sheets from an Excel File
#'
#' @description
#' Reads all sheets of a given Excel file into a named list of data frames or tibbles.
#' Each element of the list corresponds to one sheet, and is named after the sheet name.
#'
#' @param filename Character string. Path to the `.xlsx` Excel file to be read.
#'
#' @return A named list of tibbles (or data frames), one per sheet in the Excel file.
#'
#' @details
#' This function uses the `readxl` package to load each sheet in the Excel workbook.
#' It is useful when you want to batch process all sheets without specifying them manually.
#'
#' @examples
#' \dontrun{
#' all_sheets <- read_excel_allsheets("my_data.xlsx")
#' names(all_sheets) # Shows the names of the sheets
#' }
#'
#' @export
read_excel_allsheets <- function(filename) {
    # Identify all sheet names in the Excel file
    sheets <- readxl::excel_sheets(filename)

    # Read each sheet using read_excel and store in a list
    x <- lapply(sheets, function(sheetname) {
        readxl::read_excel(filename, sheet = sheetname)
    })

    # Assign sheet names as list element names
    names(x) <- sheets

    return(x)
}
