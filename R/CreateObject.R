CreateLTMGObject <- function(x = input_dir,min.cell = 0, min.gene = 0) {
  raw.matrix <- as.matrix(x)
  raw.matrix.filterbycell <- raw.matrix[(rowSums(raw.matrix > 0) > min.cell),]
  raw.matrix.filterbygene <- raw.matrix.filterbycell[(colSums(raw.matrix.filterbycell > 0) > min.gene),]
  LTMG_Object<- new(Class = 'LTMG', InputData = as.matrix(raw.matrix.filterbygene))
  return(LTMG_Object)
}



