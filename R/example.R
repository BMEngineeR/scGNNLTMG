# devtools::install_github("BMEngineeR/scGNNLTMG")

# x <- read.csv("d:/my_analysis/BRIC_TEST/2.Yan/Yan_expression.csv",header = T,row.names = 1,check.names = F)
# library(scGNNLTMG)
# object <- CreateLTMGObject(x)
# object <-RunLTMG(object,Gene_use = 100)
# WriteSparse(object,path = "d:/BMBL",gene.name = FALSE, cell.name = FALSE)
library(Seurat)
library(reticulate)
#py_install("umap-learn")
library(umap)
library(ggplot2)
setwd("d:/my_analysis/XuDong/visualization/zeisel/zeisel/")
orginal_matrix <- read.csv("original_expression.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
embedding_matrix <- read.csv("zeisel_z.csv",header = T,row.names = 1,check.names = F, stringsAsFactors = F)
label_benchmarker <- read.csv("Zeisel_7_label.csv",stringsAsFactors = F)
identical(label_benchmarker$Cell,colnames(orginal_matrix))
label_benchmarker <- label_benchmarker$Label
#label_benchmarker <- as.factor(label_benchmarker)
unique(label_benchmarker)

cell_label <- read.csv("zeisel_result.csv",header = T,check.names = F,stringsAsFactors = F)
cell_label[,1] <- colnames(orginal_matrix)
# embedding umap
rownames(embedding_matrix)<- colnames(orginal_matrix)
umap_import <- import(module = "umap", delay_load = TRUE)
set.seed(123)
embedding_umap <-  umap(
  d = embedding_matrix,
  n_neighbors = as.integer(30),
  n_components = as.integer(2),
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
  verbose = FALSE
)
my_embeding_matrix <- cbind(embedding_umap$layout[,1],embedding_umap$layout[,2],cell_label,label_benchmarker)
colnames(my_embeding_matrix)<- c("UMAP1","UMAP2","NAME","CLUSTER","Benchmark")
my_embeding_matrix$Benchmark <- as.factor(my_embeding_matrix$Benchmark)
p<-ggplot(my_embeding_matrix,aes(x=UMAP1,y=UMAP2,color = Benchmark))
p<- p+geom_point(size = 1)
p

# Seurat original imputed_matrix
MAT <- orginal_matrix
MAT <- MAT[rowSums(MAT)>0,colSums(MAT)>0]
  Gene_use_name <-rownames(MAT)[order(apply(MAT, 1, var),decreasing = T)[1:2000]]
  MAT <-MAT[Gene_use_name,]
pbmc<-CreateSeuratObject(MAT,project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
pbmc <- RunUMAP(pbmc, dims = 1:16,n.epochs = 100)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
orginal_label <- pbmc$RNA_snn_res.0.5
pbmc$celllabel <- as.factor(label_benchmarker)
Idents(pbmc)<-pbmc$celllabel
DimPlot(pbmc, reduction = "umap",pt.size = 1)

# Seurat  imputed_matrix
imputed_matrix <- read.csv("zeisel_recon.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
rownames(imputed_matrix)<- rownames(MAT)
colnames(imputed_matrix)<- colnames(MAT)
pbmc<-CreateSeuratObject(orginal_matrix,project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc@assays$RNA@data <-as.sparse(imputed_matrix)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
pbmc <- RunUMAP(pbmc, dims = 1:16,n.epochs = 100)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc$celllabel <- as.factor(label_benchmarker)
imputed_label <- pbmc$RNA_snn_res.0.5
Idents(pbmc)<-pbmc$celllabel
DimPlot(pbmc, reduction = "umap",pt.size = 1)

# ARI scGNN
library(MixGHD)
ARI(cell_label$`0`,label_benchmarker)
ARI(label_benchmarker,orginal_label)
ARI(label_benchmarker,imputed_label )

