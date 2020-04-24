#' Title
#'
#' @param x
#' @param file
#' @export
#' @return
#'
#' @examples
writeMM_new <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(c(
    sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
    sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),file)
  i = x@i+1
  j = findInterval(seq(x@x)-1,x@p[-1])+1
  value = x@x
  my.summary <- data.frame(i = i, j = j, value = value)
  write.table(my.summary,file =file,sep = " ",row.names = F,append = T,col.names = F)
}
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
writesparse <- function(object = NULL, path = getwd(), gene.name = TRUE, cell.name = TRUE) {
  if(!require("Matrix")){
    install.packages("Matrix")
  }
  my.sparse <- Matrix(object@OrdinalMatrix,sparse = T)
  if(all(my.sparse <=1)){
    writeMM_new(x=my.sparse,file = paste0(path,"/LTMG_sparse.mtx"))
  }else{
    writeMM(my.sparse,file = paste0(path,"/LTMG_sparse.mtx"))
  }
  if(cell.name == TRUE){
    write.table(data.frame(Barcode =colnames(my.sparse)) ,file= paste0(path,"/barcode.txt"),sep = "\t",quote = F,row.names = F)
  }
  if(gene.name == TRUE){
    write.table(data.frame(Gene = rownames(my.sparse)),file = paste0(path,"/Gene.txt"),sep = "\t",quote = F,row.names = F)
  }
}

#' @export
#' @rdname WriteSparse
setMethod("WriteSparse","LTMG", writesparse)
