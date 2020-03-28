#
# setwd("d:/my_analysis/BRIC_TEST/2.Yan/")
# x <- read.csv("d:/my_analysis/BRIC_TEST/2.Yan/Yan_expression.csv",header = T,row.names = 1,check.names = F)
# library(scGNNLTMG)
# object <- CreateLTMGObject(x)
# object <-RunLTMG(object,Gene_use = 100)
#
# my.sparse <- Matrix::Matrix(object@OrdinalMatrix,sparse = T)
# writeMM(my.sparse,file = "LTMG_sparse.mtx")
# write.table(data.frame(Barcode =colnames(my.sparse)) ,file= "barcode.txt",sep = "\t",quote = F,row.names = F)
# write.table(data.frame(Gene = rownames(my.sparse)),file = "Gene.txt",sep = "\t",quote = F,row.names = F)
