
#' Title
#'
#' @param object
#'
#' @return
#' @export
#' @name WriteSparse
#' @examples
#' @importFrom Matrix Matirx writeMM
writesparse <- function(object = NULL) {
  my.sparse <- Matrix(object@OrdinalMatrix,sparse = T)
  writeMM(my.sparse,file = "LTMG_sparse.mtx")
  write.table(data.frame(Barcode =colnames(my.sparse)) ,file= "barcode.txt",sep = "\t",quote = F,row.names = F)
  write.table(data.frame(Gene = rownames(my.sparse)),file = "Gene.txt",sep = "\t",quote = F,row.names = F)

}

#' @export
#' @rdname WriteSparse
setMethod("WriteSparse","LTMG", writesparse)
