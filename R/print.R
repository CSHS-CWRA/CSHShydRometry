# The print method shared by every rating-curve class.

#' Print a rating curve fit
#'
#' Every fit carries `"rating_curve"` as its second class, so this one method
#' covers all of them.
#'
#' @param x A rating curve fit.
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the printed output.
#' @export
print.rating_curve <- function(x, ...) {
  cl <- class(x)
  cl <- cl[cl != "rating_curve"]
  cl <- paste(cl, collapse = ", ")
  cat("Rating curve model.\n- Method:", cl, "\n")
}
