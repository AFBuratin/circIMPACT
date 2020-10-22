#'
#'
#' @param countmatrix
#'
#' @return list
#'
#' @examples
#' dataset <- matrix(sample(c(NA, 1:5), 25, replace = TRUE), 5)
#' df <- as.data.frame(dataset)
#' package::nazero(df)
#'
#' @export

# The ‘export’ field will contain the function names that you want the end user to access.
# You can also hide certain functions from your code that you don’t want to offer by not mentioning them here.

nazero <- function(countmatrix)
{
  countmatrix[is.na(countmatrix)] <- 0
  return(countmatrix)
}
