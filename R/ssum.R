#' A helping function to add array elements of a list
#'
#' @param .list a list of three-dimensional arrays to be cumulatively summed
#' @keywords preparing data
#' @export 
#' @examples
#' ssum()
 
ssum=function(.list) {
  if(is.array(.list)) .list=array.to.list(.list)
  return(Reduce("+", .list))
}
