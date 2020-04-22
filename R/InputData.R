
#' Title
#'
#' @param object
#' @param path output path
#' @param gene.name whether output gene name for sparse matrix
#' @param cell.name whether output cell name for sparse matrix
#'
#' @return
#' @name WriteSparse
#' @examples
#' @importFrom  Matrix Matrix writeMM
writesparse <- function(object = NULL, path = "./", gene.name = TRUE, cell.name = TRUE) {
  if(!require("Matrix")){
    install.packages("Matrix")
  }
  my.sparse <- Matrix(object@OrdinalMatrix,sparse = T)
  writeMM(my.sparse,file = paste0(path,"/LTMG_sparse.mtx"))
  if(cell.name == TRUE){
    write.table(data.frame(Barcode =colnames(my.sparse)) ,file= paste0(path,"barcode.txt"),sep = "\t",quote = F,row.names = F)
  }
  if(gene.name == TRUE){
    write.table(data.frame(Gene = rownames(my.sparse)),file = paste0(path,"Gene.txt"),sep = "\t",quote = F,row.names = F)
  }
}

#' @export
#' @rdname WriteSparse
setMethod("WriteSparse","LTMG", writesparse)
