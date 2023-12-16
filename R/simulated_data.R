#' Simulated PDF data.
#'
#' Data simulated using the the dfnWorks suite.  
#'
#' `pdf_data` is a data frame with 1,000 rows and 5 columns.
#' `x_grid`; is a timestamp of each of `full_data`'s 1,000 rows.
#'
#' @format `pdf_data` is an n x p matrix, where n is the number of PDFs and p matches the length of x_grid.  x_grid contains the points at which the PDFs are evaluated (assumed equally spaced apart).
#' @examples
#' pdf_data
#' x_grid
"pdf_data"

#' @rdname pdf_data
#' @format NULL
"x_grid"
