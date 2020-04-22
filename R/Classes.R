setClass("LTMG",slots = c(
  InputData = "ANY",
  OrdinalMatrix = "matrix"
)
)

setClass("Vis",slots = c(
  OriginalMatrix = "ANY",
  Cluster = "ANY",
  Embedding = "ANY",
  Reduction = "ANY",
  Network= "ANY",
  ImputatedData = "ANY",
  MarkerGene = "ANY",
  GeneInfo = "ANY",
  CellInfo ="ANY",
  activated.idents = "ANY",
  tmp.seurat = "ANY"
))
