# Copilot instructions (Markerscreen-Tcell)

## R style & packaging

- Use **roxygen2 tags** explicitly:
  - Include `@title` and `@description` for exported functions.
  - Add `@details` for scientific/background context and assumptions.
  - Add educational `@examples` and include short `# Expected:` comments when helpful.
- Use vertical space economically. I.e. write if-else and or for loops in a
  single line, if it can be done in less than ~80 characters.
- Try to fit roxygen parameters in a single line with less than ~80 characters
  in total. If not possible, use multiple lines with the description starting on
  the next line. Examples:
- Function should have less than 50 lines. Ideal is a maximum of 10-25 lines.
  Functions shorter than 10 lines must not be used unless the function is called
  at least twice in the codebase. Functions between 10 and 50 lines can and should
  also be used for logical code separation, even if only called once. This makes
  the code more readable and testable.

  ```R
  ## Good:
  #' @param x A numeric vector of input values.

  ## Good:
  #' @param long_param_name
  #' A character string that specifies the name of the output file to write.
  #' Environment variables in the form ${VAR} will be expanded.
  #'

  ## Bad:
  #' @param long_param_name A character string that specifies the name of the
  #'   output file to write. Environment variables in the form ${VAR} will
  #'   be expanded.

  ## Good
  #' @return
  #' A data.frame with columns file_path, round, plate, marker, well,
  #' harmonized, and outfile.
  #
  ## Bad
  #' @return A data.frame with columns file_path, round, plate, marker,
  #'   well,harmonized, and outfile.
  ```


## Dependencies

- Prefer **base R** (and the packages that ship with R by default) over tidyverse.
- Avoid the pipe (`|>`, `%>%`) in new code when practical.
- Avoid adding new tidyverse dependencies; if parsing/reshaping can be done in base R, do so.
- Cytometry-specific packages (e.g., `flowCore`, `flowWorkspace`, `cytoKal`, `cytobench`) are OK to use when required.

## Design expectations

- Functions should have explicit inputs/outputs; avoid reliance on global variables.
- Validate inputs early with clear error messages.
- Write outputs only under `res/`.

## CI/CD expectations

- All new code should include unit tests where practical.
- All outputs must be written under `res/`.
- Documentation should be updated for any new exported functions.
- The main pipeline is orchestrated via `main.R` using `devtools::load_all()` and function calls (no script sourcing).

For more details on style and expectations, see `.github/copilot-instructions.md`.