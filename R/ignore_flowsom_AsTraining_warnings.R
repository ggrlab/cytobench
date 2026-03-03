ignore_flowsom_AsTraining_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)

      if (msg %in% c(
        "Data will be compensated as was done while computing the FlowSOM object.",
        "Data will be transformed as was done while computing the FlowSOM object.",
        "Data will be scaled as was done while computing the FlowSOM object."
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
