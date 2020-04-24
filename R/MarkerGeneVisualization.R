
#' @import ggplot2
#' Title
#'
#' @param object
#' @param dim
#'
#' @return
#' @export
#'
#' @examples
FindMarkerGene <- function(object = NULL,min.pct = 0.25, logfc.threshold = 0.25){
  tmp.seurat <- Seurat::CreateSeuratObject(t(object@ImputatedData))
  tmp.seurat <- Seurat::FindVariableFeatures(tmp.seurat, selection.method = "vst", nfeatures = length(object@GeneInfo))
  tmp.seurat$cell_label <- object@Cluster$cluster
  Seurat::Idents(tmp.seurat) <-  tmp.seurat$cell_label
  object@MarkerGene  <- Seurat::FindAllMarkers(tmp.seurat, only.pos = TRUE, test.use = "wilcox", min.pct = min.pct, logfc.threshold = logfc.threshold)
  object@tmp.seurat <- tmp.seurat
  return(object)
}



PlotHeatmap <- function(object = NULL,top.gene = 10,p.adj = 0.05){
  marker.list <- object@MarkerGene
  marker.list <- marker.list[marker.list$p_val_adj < p.adj,]
  marker.cluster.index <- as.data.frame(cbind(index =  1:length(unique(marker.list$cluster)), cluster= unique(as.character(marker.list$cluster))),
                                        stringsAsFactors = F)
  sub.marker.list<- c()
  for (i in 1:length(unique(marker.list$cluster))){
    tmp.cluster <- marker.cluster.index$cluster[marker.cluster.index$index == i]
    tmp.marker.list <- marker.list[marker.list$cluster == tmp.cluster,]
    tmp.marker.list <- tmp.marker.list[1:top.gene,]
    sub.marker.list <- rbind(sub.marker.list,tmp.marker.list)
  }
  my.cluster<- object@Cluster
  my.cluster <- my.cluster[order(my.cluster$cluster),]
  if(!is.factor(my.cluster$cluster)){
    col.data <- as.factor(my.cluster$cluster)
  } else{col.data <- my.cluster$cluster}
  heatmap.matrix <- object@OriginalMatrix[sub.marker.list$gene,as.character(my.cluster$cell)]
  heatmap.matrix <- na.omit(heatmap.matrix)
  if(!require(colorspace)){
    install.packages("colorspace")
  }
  colSide <- qualitative_hcl(length(unique(col.data)),palette = "Dynamic")[col.data]
  colMain <- colorRampPalette(c("#2D80BD","#DCEAF5","#B0390C"))(100)
  heatmap(as.matrix(heatmap.matrix), Colv = NA, Rowv = NA, scale = "column",ColSideColors =colSide ,col = colMain,labCol = "0")
  legend("topright",legend=c(unique(as.character(my.cluster$cluster))), bg="transparent",pch = rep(20,length(unique(col.data))),
         col=qualitative_hcl(length(unique(col.data)),palette = "Dynamic")[unique(my.cluster$cluster)], cex=0.8, box.lty=0)
}

PlotViolin <- function(object = NULL,feature.name = NULL, pt.show =FALSE){
  if(!any(grepl(paste0("^",feature.name,"$"),rownames(object@OriginalMatrix),ignore.case = T))){
    stop("Does not find gene ", feature.name, ". Maybe: ", paste0(grep(feature.name,rownames(object@OriginalMatrix),ignore.case = T,value = T),collapse = ","))
  }
  index <- grep(paste0("^",feature.name,"$"),rownames(object@OriginalMatrix),ignore.case = T)
  if(!identical(colnames(object@OriginalMatrix),object@CellInfo)){
    object@OriginalMatrix <- object@OriginalMatrix[,match(object@CellInfo,colnames(object@OriginalMatrix))]
  }
  violin.matrix <- as.data.frame(cbind(ID= object@CellInfo, VALUE= as.numeric(object@OriginalMatrix[index,]), CLUSTER = as.character(object@Cluster$cluster)))
  rownames(violin.matrix) <- violin.matrix$ID
  if(!require(colorspace)){
    install.packages("colorspace")
  }
  colSide <- qualitative_hcl(length(unique(violin.matrix$CLUSTER)),palette = "Dynamic")[violin.matrix$CLUSTER]
  color.idx <- cbind(CLUSTER= unique(as.character(violin.matrix$CLUSTER)),values=unique(colSide))
  color.idx <- color.idx[order(color.idx[,1]),]
  p <- ggplot(violin.matrix,aes(x = as.factor(CLUSTER), y = as.numeric(as.character(VALUE))))
  if(pt.show == TRUE){
  p  <- p+ geom_violin()+ geom_jitter(height = 0, width = 0.05) +theme_classic()
  }
  p  <- p+ geom_violin(aes(fill=factor(CLUSTER)))+theme_classic() + scale_fill_manual(breaks = color.idx[,1],values=color.idx[,2])
  p <- p+ labs(title =feature.name ,x = "Cluster", y = "Value",fill = "cluster")
  print(p)
}

PlotScatter <- function(object = NULL, feature.name = NULL){
  if(length(feature.name)!=2){
    stop("Only accept two genes, please provide the correct gene number.")
  }
  first.gene <- feature.name[1]
  second.gene <- feature.name[2]
  if(!any(grepl(paste0("^",first.gene,"$"),rownames(object@OriginalMatrix),ignore.case = T))){
    stop("Does not find gene ", first.gene, ". Maybe: ", paste0(grep(first.gene,rownames(object@OriginalMatrix),ignore.case = T,value = T),collapse = ","))
  }
  if(!any(grepl(paste0("^",second.gene,"$"),rownames(object@OriginalMatrix),ignore.case = T))){
    stop("Does not find gene ", second.gene, ". Maybe: ", paste0(grep(second.gene,rownames(object@OriginalMatrix),ignore.case = T,value = T),collapse = ","))
  }
  index.1 <- grep(paste0("^",first.gene,"$"),rownames(object@OriginalMatrix),ignore.case = T)
  index.2 <- grep(paste0("^",second.gene,"$"),rownames(object@OriginalMatrix),ignore.case = T)
  if(!identical(colnames(object@OriginalMatrix),object@CellInfo)){
    object@OriginalMatrix <- object@OriginalMatrix[,match(object@CellInfo,colnames(object@OriginalMatrix))]
  }
  scatter.matrix <- as.data.frame(cbind(ID= object@CellInfo,
                                       VALUE1= as.numeric(object@OriginalMatrix[index.1,]),
                                       VALUE2= as.numeric(object@OriginalMatrix[index.2,]),
                                       CLUSTER = object@Cluster$cluster))
  rownames(scatter.matrix)  <- scatter.matrix$ID
  R.value <-round(cor(as.numeric(as.character(scatter.matrix$VALUE1)),as.numeric(as.character(scatter.matrix$VALUE2)),method ="pearson" ),3)
  p <- ggplot(scatter.matrix,aes(x=as.numeric(as.character(VALUE1)),y=as.numeric(as.character(VALUE2))))
  p <- p + geom_point() +labs(title = paste0("Pearson R = ",R.value),x = first.gene,y = second.gene)
  p <- p + geom_smooth(method=lm, fullrange=TRUE) + theme_classic()
  p
}

##








