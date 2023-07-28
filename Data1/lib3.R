#Load libraries
library(Seurat)
library(patchwork)
library(clustree)
library(DoubletFinder)
library(dplyr)    # alternatively, this also loads %>%
# Load the data 
lib3.data <- Read10X("/home/yanis/Bureau/Finale/Data1/larval3/filtered_feature_bc_matrix/")
lib3 <-CreateSeuratObject(counts =lib3.data, min.cells = 3, min.features = 200,project="Lib3")
VlnPlot(lib3, features = c("nFeature_RNA", "nCount_RNA") , ncol = 2)
lib3 <- subset(lib3 ,subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
lib3 <- NormalizeData(object = lib3, verbose = FALSE)
lib3 <- FindVariableFeatures(object = lib3, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
lib3 <- ScaleData(object = lib3)
lib3 <- RunPCA(object = lib3 )
ElbowPlot(lib3 , ndims = 50)
lib3
lib3 <- FindNeighbors(object = lib3, dims = 1:50)
lib3 <- FindClusters(object = lib3)
lib3 <- RunUMAP(object = lib3, dims = 1:50)




## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lib3 <- paramSweep_v3(lib3, PCs = 1:50, sct = FALSE)
sweep.res.list_lib3
sweep.stats_lib3 <- summarizeSweep(sweep.res.list_lib3, GT = FALSE)
sweep.stats_lib3
bcmvn_lib3 <- find.pK(sweep.stats_lib3)

ggplot(bcmvn_lib3, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_lib3 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))
pK  # 0.02

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib3@meta.data$seurat_clusters
annotations
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
homotypic.prop
nExp_poi <- round(0.13*nrow(lib3@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
lib3 <- doubletFinder_v3(lib3,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
lib3@meta.data
pdf("visualuzation_D_Slib3.pdf")
DimPlot(lib3, reduction = 'umap', group.by = "DF.classifications_0.25_0.02_1666")
dev.off()


pdf("dpn_distrubutionlib3.pdf")
FeaturePlot(lib3 , features = "LOC6627282")
dev.off()

# number of singlets and doublets
table(lib3@meta.data$DF.classifications_0.25_0.02_1666)
#Doublet Singlet 
#1666   12964

object_NoDoublets <- subset(lib3, subset =  DF.classifications_0.25_0.02_1666== "Singlet")
object_NoDoublets
saveRDS(object_NoDoublets , "lib3_NoDoublet.rds")


object_NoDoublets <- readRDS("lib3_NoDoublet.rds")
object_NoDoublets









