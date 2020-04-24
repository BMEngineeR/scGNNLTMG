PlotNetwork <- function(object = NULL,seed = 123){
  if(!require(GGally)){
    message("installing required package: GGally")
    install.packages("GGally")
  }
  if(!require(devtools)){
    message("installing required package: devtools")
    install.packages("devtools")
  }
  if(!require(ggnet)){
    message("installing required package: ggnet")
    devtools::install_github("briatte/ggnet")
  }
  if(!require(network)){
    message("installing required package: network")
    install.packages("network")
  }
  if(!require(sna)){
    message("installing required package: sna")
    install.packages("sna")
  }
  if(!require(Matrix)){
    message("installing required package: Matrix")
    install.packages("Matrix")
  }
  my.network.data<- object@Network
  my.network.matrix <- sparseMatrix(i = as.numeric(as.character(my.network.data$node1)) + 1,
                                    j =  as.numeric(as.character(my.network.data$node2)) +1,
                                    x =  as.numeric(as.character(my.network.data$edge)),
                                    dims = c(max(as.numeric(as.character(my.network.data$node1)))+1, max(as.numeric(as.character(my.network.data$node2)))+1))
  my.network.matrix <- as.matrix( my.network.matrix )
  rownames(my.network.matrix) <- unique(my.network.data$cell_info)
  colnames(my.network.matrix) <- unique(my.network.data$cell_info)
  net = network(as.matrix(my.network.matrix), directed = FALSE)
  if(!require(colorspace)){
    install.packages("colorspace")
  }
  col.data <- object@Cluster$cluster
  colSide <- qualitative_hcl(length(unique(col.data)),palette = "Dynamic")[col.data]
  color.idx <- cbind(CLUSTER= unique(as.character(col.data)),values=unique(colSide))
  color.idx <- color.idx[order(color.idx[,1]),]
  set.seed(seed)
  ggnet2(net,node.size = 2,color = colSide)
}

PlotNetworkOneCluster <- function(object = NULL,cluster.idx = "0",seed = 123){
  if(!require(GGally)){
    message("installing required package: GGally")
    install.packages("GGally")
  }
  if(!require(devtools)){
    message("installing required package: devtools")
    install.packages("devtools")
  }
  if(!require(ggnet)){
    message("installing required package: ggnet")
    devtools::install_github("briatte/ggnet")
  }
  if(!require(network)){
    message("installing required package: network")
    install.packages("network")
  }
  if(!require(sna)){
    message("installing required package: sna")
    install.packages("sna")
  }
  if(!require(Matrix)){
    message("installing required package: Matrix")
    install.packages("Matrix")
  }
  my.network.data<- object@Network
  my.network.matrix <- sparseMatrix(i = as.numeric(as.character(my.network.data$node1)) + 1,
                                    j =  as.numeric(as.character(my.network.data$node2)) +1,
                                    x =  as.numeric(as.character(my.network.data$edge)),
                                    dims = c(max(as.numeric(as.character(my.network.data$node1)))+1, max(as.numeric(as.character(my.network.data$node2)))+1))
  my.network.matrix <- as.matrix( my.network.matrix )
  rownames(my.network.matrix) <- unique(my.network.data$cell_info)
  colnames(my.network.matrix) <- unique(my.network.data$cell_info)
  if(!as.character(cluster.idx) %in% as.character(my.network.data$node1_cluster)){
    stop("cannot find cluster " ,cluster.idx,
         " in this cluster information.",
         " please choose one of ",
         paste0(unique(as.character(my.network.data$node1_cluster)),collapse = ", ")
         )
  }
  # select cell according to idx
  select.cell <- object@Cluster$cell[as.character(object@Cluster$cluster) %in% as.character(cluster.idx)]
  my.network.submatrix <- my.network.matrix[ select.cell, select.cell]
  net.sub = network(as.matrix(my.network.submatrix), directed = FALSE)
  if(!require(colorspace)){
    install.packages("colorspace")
  }
  col.data <- object@Cluster$cluster
  colSide <- qualitative_hcl(length(unique(col.data)),palette = "Dynamic")[col.data]
  color.idx <- cbind(CLUSTER= unique(as.character(col.data)),values=unique(colSide))
  color.idx <- color.idx[order(color.idx[,1]),]
  color.idx <- as.data.frame(color.idx)
  # select cluster specific color
  color.select <- as.character(color.idx$values)[as.character(color.idx$CLUSTER) %in% as.character(cluster.idx)]
  set.seed(seed)
  p<- ggnet2(net.sub,node.size = 2,color = rep(color.select,nrow(my.network.submatrix)))
  p <- p+ labs(title = paste0("Cell similarity in cluster ",as.character(cluster.idx)))
  print(p)
}


PlotNetwork_sampling <- function(object = NULL,seed = 123){
  my.network.data<- object@Network
  my.network.matrix <- sparseMatrix(i = as.numeric(as.character(my.network.data$node1)) + 1,
                                    j =  as.numeric(as.character(my.network.data$node2)) +1,
                                    x =  as.numeric(as.character(my.network.data$edge)),
                                    dims = c(max(as.numeric(as.character(my.network.data$node1)))+1, max(as.numeric(as.character(my.network.data$node2)))+1))
  my.network.matrix <- as.matrix( my.network.matrix )
  rownames(my.network.matrix) <- unique(my.network.data$cell_info)
  colnames(my.network.matrix) <- unique(my.network.data$cell_info)
  my.subsample<- c()
  for (i in 1:length(unique(as.character(my.network.data$node1_cluster)))){
    tmp.cell.cluster <- unique(as.character(my.network.data$node1_cluster))[i]
    tmp.cell.name <- object@Cluster$cell[object@Cluster$cluster %in%  tmp.cell.cluster]
    set.seed(123)
    tmp.sample.name <- tmp.cell.name[sample(1:length(tmp.cell.name),50)]
    my.subsample <- c(my.subsample,tmp.sample.name)
}
  my.sub.matrix <- as.matrix(my.network.matrix)[my.subsample,my.subsample]
  net = network(my.sub.matrix, directed = FALSE)
  if(!require(colorspace)){
    install.packages("colorspace")
  }
  col.data <- object@Cluster$cluster[object@Cluster$cell %in% rownames(my.sub.matrix)]
  colSide <- qualitative_hcl(length(unique(col.data)),palette = "Dynamic")[col.data]
  color.idx <- cbind(CLUSTER= unique(as.character(col.data)),values=unique(colSide))
  color.idx <- color.idx[order(color.idx[,1]),]
  set.seed(seed)
  ggnet2(net,node.size = 2,color = colSide)
}
