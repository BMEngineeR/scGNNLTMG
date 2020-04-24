
RunUmap <- function(object = NULL,
                     seed = 123,
                     n_neighbors = 30L,
                     n_components = 2L,
                     metric = "cosine",
                     n_epochs = 100,
                     learning_rate = 1.0,
                     min_dist = 0.3,
                     spread = 1.0,
                     set_op_mix_ratio = 1.0,
                     local_connectivity = 1L ,
                     repulsion_strength = 1,
                     negative_sample_rate = 5,
                     a = NA,
                     b = NA,
                     fast_sgd = FALSE,
                     verbose = FALSE){
  # check required package
  if(!require(umap)){
    message("installing required package: umap")
    install.packages("umap")
  }
  if(!require(ggplot2)){
    message("installing required package: ggplot2")
    install.packages("ggplot2")
  }
  embedding_matrix <- object@Embedding
  set.seed(seed)
  embedding_umap <-  umap(
    d = embedding_matrix,
    n_neighbors = as.integer(n_neighbors),
    n_components = as.integer(n_components),
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    min_dist = min_dist,
    spread = spread,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity ,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    fast_sgd = fast_sgd,
    verbose =verbose
  )
  object@Reduction <- embedding_umap
  return(object)
  }
PlotUmap <- function(obejct){
  my_embeding_matrix<- object@Reduction
  label_cell <- object@Cluster
  my_embeding_matrix <- cbind(my_embeding_matrix$layout[,1],my_embeding_matrix$layout[,2],label_cell)
  colnames(my_embeding_matrix)[c(1,2)]<- c("UMAP1","UMAP2")
  my_embeding_matrix$cluster <- as.factor(my_embeding_matrix$cluster)
  if(!require(colorspace)){
    install.packages("colorspace")
  }
  colSide <- qualitative_hcl(length(unique(my_embeding_matrix$cluster)),palette = "Dynamic")[my_embeding_matrix$cluster]
  color.idx <- cbind(cluster= unique(as.character(my_embeding_matrix$cluster)),values=unique(colSide))
  color.idx <- color.idx[order(color.idx[,1]),]
  p<-ggplot(my_embeding_matrix,aes(x=UMAP1,y=UMAP2,color = cluster ))
  p<- p+geom_point(size = 1)+theme_classic()
  p <- p + scale_color_manual(breaks = color.idx[,1],values=color.idx[,2])
  print(p)
}

PlotGenes <- function(object, feature.name = NULL){
  my_embeding_matrix<- object@Reduction
  my_embeding_matrix <- cbind(my_embeding_matrix$layout[,1],my_embeding_matrix$layout[,2])
  colnames(my_embeding_matrix)[c(1,2)]<- c("UMAP1","UMAP2")
  my.impute.matrix <- t(object@ImputatedData)
  if(!any(grepl(paste0("^",feature.name,"$"),rownames(my.impute.matrix),ignore.case = T))){
    stop("Does not find gene in imputate matrix", feature.name, ". Maybe: ", paste0(grep(feature.name,rownames(my.impute.matrix),ignore.case = T,value = T),collapse = ","))
  }
  index <- grep(paste0("^",feature.name,"$"),rownames(my.impute.matrix),ignore.case = T)
  my_embeding_matrix <- as.data.frame(my_embeding_matrix)
  my_embeding_matrix$VALUE <- my.impute.matrix[index,rownames(my_embeding_matrix)]
  my_embeding_matrix <- as.data.frame(my_embeding_matrix)
  p<-ggplot(my_embeding_matrix,aes(x=UMAP1,y=UMAP2,color = VALUE ))
  p<- p+geom_point(size = 1)+theme_classic()
  p <- p + scale_color_gradient(low = "#BFC2C8" ,high = "#FF0F04")
  p <- p + labs(title = paste0(feature.name," imputed expression") ,x = "UMAP1", y = "UMAP2",fill = "Value")
  print(p)

}

