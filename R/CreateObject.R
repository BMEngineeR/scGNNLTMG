CreateLTMGObject <- function(x = input_dir,min.cell = 0, min.gene = 0) {
  raw.matrix <- as.matrix(x)
  raw.matrix.filterbycell <- raw.matrix[(rowSums(raw.matrix > 0) > min.cell),]
  raw.matrix.filterbygene <- raw.matrix.filterbycell[,(colSums(raw.matrix.filterbycell > 0) > min.gene)]
  message("Creating BRIC object. \n",
          "The original input file contains ", dim(raw.matrix)[2], " cells and ", dim(raw.matrix)[1], " genes \n",
          "Removed ", dim(raw.matrix)[1] - dim(raw.matrix.filterbycell)[1], " genes that total expression value is equal or less than ", min.cell, "\n",
          "Removed ", dim(raw.matrix.filterbycell)[2] - dim(raw.matrix.filterbygene)[2], " cells that number of expressed gene is equal or less than ", min.gene
  )
  LTMG_Object<- new(Class = 'LTMG', InputData =  raw.matrix.filterbygene)
  return(LTMG_Object)
}



