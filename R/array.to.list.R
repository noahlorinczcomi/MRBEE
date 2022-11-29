#' A helping function to convert a 3d array into a list
#'
#' @param .array an array of dimensions axbxc to be converted into a list of length c with elements axb matrices
#' @keywords preparing data
#' @export 
#' @examples
#' array.to.list()
 
array.to.list=function(.array) lapply(seq(dim(.array)[3]), function(x) .array[ , , x])
