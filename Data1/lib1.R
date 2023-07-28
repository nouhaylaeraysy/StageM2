#Load libraries
library(Seurat)
library(patchwork)
library(clustree)
library(DoubletFinder)
library(dplyr)    # alternatively, this also loads %>%

# Load the data 
lib1.data <- Read10X("/home/yanis//Bureau/Finale/Data1/larval1/filtered_feature_bc_matrix/")
lib1 <-CreateSeuratObject(counts =lib1.data, min.cells = 3, min.features = 200,project="Lib1")
lib1
VlnPlot(lib1, features = c("nFeature_RNA", "nCount_RNA") , ncol = 2)
lib1 <- subset(lib1 ,subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
lib1 <- NormalizeData(object = lib1, verbose = FALSE)
lib1 <- FindVariableFeatures(object = lib1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
lib1 <- ScaleData(object = lib1)
lib1 <- RunPCA(object = lib1 , npcs = 100 )
ElbowPlot(lib1 , ndims = 100)
lib1
lib1 <- FindNeighbors(object = lib1, dims = 1:50)
lib1 <- FindClusters(object = lib1)
lib1 <- RunUMAP(object = lib1, dims = 1:50)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lib1 <- paramSweep_v3(lib1, PCs = 1:50, sct = FALSE)
sweep.res.list_lib1
sweep.stats_lib1 <- summarizeSweep(sweep.res.list_lib1, GT = FALSE)
sweep.stats_lib1
bcmvn_lib1 <- find.pK(sweep.stats_lib1)

ggplot(bcmvn_lib1, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_lib1 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))
pK  # 0.3

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib1@meta.data$seurat_clusters
annotations
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
homotypic.prop
nExp_poi <- round(0.13*nrow(lib1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
lib1 <- doubletFinder_v3(lib1,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)
lib1@meta.data

# visualize doublets
pdf("visualuzation_D_S.pdf")
DimPlot(lib1, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_1540")
dev.off()
pdf("dpn_distrubutionlib1.pdf")
FeaturePlot(lib1 , features = "LOC6627282")
dev.off()
# number of singlets and doublets
table(lib1@meta.data$DF.classifications_0.25_0.3_1540)
#Doublet Singlet 
#1540   11719

object_NoDoublets <- subset(lib1, subset =  DF.classifications_0.25_0.3_1540== "Singlet")
object_NoDoublets
saveRDS(object_NoDoublets , "lib1_NoDoublet.rds")

object_NoDoublets <- readRDS("lib1_NoDoublet.rds")
object_NoDoublets










