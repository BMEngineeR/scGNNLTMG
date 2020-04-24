CreateLTMGObject <- function(x = NULL,min.cell = 0, min.gene = 0) {
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

#' Title
#'
#' @param input_path input directory path which contain three major output and original matrix.
#' @param original_matrix input original expression matrix
#' @param min.cell
#' @param min.gene
#'
#' @return
#' @export
#'
#' @examples
CreateVisObject <- function(original_matrix = NULL, input_path = NULL) {
  input_path_result <- paste0(input_path,"/results")

  # import used gene
  gene_file <-list.files(path = input_path, pattern = "Gene.txt")
  gene_input <- file.path(input_path,gene_file)
  gene_used <- read.table(gene_input,stringsAsFactors = F,header = T)
  # import cell name
  cell_file <- list.files(path = input_path, pattern = "barcode.txt")
  cell_input <- file.path(input_path,cell_file)
  cell_used <- read.table(cell_input,stringsAsFactors = F,header = T)

  # import embedding
  embedding_file <-  list.files(path =  input_path_result,pattern = "embedding.csv")
  embedding_input <-file.path( input_path_result,embedding_file)
  embedding_matrix <- read.csv(embedding_input,header = F,check.names = F, stringsAsFactors = F)
  rownames(embedding_matrix) <- cell_used$Barcode

  # import cell label
  label_file <- list.files(path =  input_path_result,pattern = "results.txt")
  label_input <- file.path( input_path_result,label_file)
  label_cell <- read.table(label_input,stringsAsFactors = F,header = F,check.names = F)
  label_cell <- as.data.frame(cbind(ID = rownames(embedding_matrix),CLUSTER=label_cell$V1))
  colnames(label_cell) <- c("cell","cluster")
  label_cell$cell <- rownames(embedding_matrix)
  # import network
  network_file <- list.files(path = input_path_result,pattern = "graph.csv")
  network_input <- file.path(input_path_result,network_file)
  network_matrix <- read.csv(network_input,stringsAsFactors = F,header = F)
  cluster.label <- read.table(label_input,stringsAsFactors = F,header = F,check.names = F)
  cluster.label$V2 <- 0:(nrow(label_cell)-1)
  cluster.label$V3 <- rownames(embedding_matrix)
  node1 <- network_matrix[,1]
  node2 <- network_matrix[,2]
  edge <- network_matrix[,3]
  node1_cluster <- rep(NA, nrow(network_matrix))
  node2_cluster <- rep(NA, nrow(network_matrix))
  cell_info <- rep(NA,nrow(network_matrix))
  for (i in 1:nrow(cluster.label)){
    node1_cluster <- ifelse(node1==cluster.label$V2[i],cluster.label$V1[i],node1_cluster)
    node2_cluster <- ifelse(node2==cluster.label$V2[i],cluster.label$V1[i],node2_cluster)
    cell_info <- ifelse(node1==cluster.label$V2[i],cluster.label$V3[i],cell_info)
  }
  network_matrix_data <- cbind(node1 = node1, node2 = node2, edge = edge,  cell_info=  cell_info,node1_cluster = node1_cluster,node2_cluster= node2_cluster)
  # import imputated matrix
  impute_file <- list.files(path =  input_path_result,pattern = "recon.csv")
  impute_input <- file.path(input_path_result,impute_file)
  impute_file <- read.csv(impute_input,header = F,stringsAsFactors = F,check.names = F)
  colnames(impute_file) <- gene_used$Gene
  rownames(impute_file) <- cell_used$Barcode
  # import orginal
  original_matrix <- original_matrix
  Vis_Object<- new(Class = 'Vis',
                    OriginalMatrix =  original_matrix,
                    Cluster = label_cell,
                    Network = as.data.frame(network_matrix_data),
                    Embedding = embedding_matrix,
                    ImputatedData =impute_file,
                    CellInfo = cell_used$Barcode,
                    GeneInfo = gene_used$Gene)
  return(Vis_Object)
}

SwitchLabel <- function(object, switch_label = NULL ){
  my.cluster<-object@Cluster
  switch_label <- switch_label
  switch.cluster <- rep(NA,nrow(object@Cluster))
  for(i in 1:nrow(switch_label)){
    tmp.original.label <- names(switch_label[i,])
    switch.cluster <- ifelse(as.character(my.cluster$cluster)== as.character(tmp.original.label),
                             as.character(switch_label[i]),switch.cluster)
  }
  my.cluster.new<-cbind(cell=my.cluster[,1],cluster = as.character(switch.cluster),cluster_backup = as.character(my.cluster[,2]))
  object@Cluster <- as.data.frame(my.cluster.new)
  return(object)
}
