
#' Title
#'
#' @param object
#'
#' @return
#' @name WriteSparse
#' @examples
#' @importFrom  Matrix Matrix writeMM
writesparse <- function(object = NULL, path = "./") {
  if(!require("Matrix")){
    install.packages("Matrix")
  }
  my.sparse <- Matrix(object@OrdinalMatrix,sparse = T)
  writeMM(my.sparse,file = paste0(path,"/LTMG_sparse.mtx"))
  write.table(data.frame(Barcode =colnames(my.sparse)) ,file= paste0(path,"barcode.txt"),sep = "\t",quote = F,row.names = F)
  write.table(data.frame(Gene = rownames(my.sparse)),file = paste0(path,"Gene.txt"),sep = "\t",quote = F,row.names = F)

}

#' @export
#' @rdname WriteSparse
setMethod("WriteSparse","LTMG", writesparse)
