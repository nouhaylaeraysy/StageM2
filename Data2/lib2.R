#Load libraries 
library(Seurat)
library(patchwork)
library(clustree)
library(DoubletFinder)
library(dplyr)    # alternatively, this also loads %>%
# Load the data 
lib2.data <- Read10X("/home/yanis/Bureau/Finale/filtered_feature_bc_matrix_larva2_virilis/")
lib2 <-CreateSeuratObject(counts =lib2.data, min.cells = 3, min.features = 200,project="Lib2")
lib2
VlnPlot(lib2, features = c("nFeature_RNA", "nCount_RNA") , ncol = 2)
lib2 <- subset(lib2 ,subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
lib2 <- NormalizeData(object = lib2, verbose = FALSE)
lib2 <- FindVariableFeatures(object = lib2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
lib2 <- ScaleData(object = lib2)
lib2 <- RunPCA(object = lib2 , npcs = 100 )
lib2
ElbowPlot(lib2 , ndims = 100)
lib2 <- FindNeighbors(object = lib2, dims = 1:50)
lib2 <- FindClusters(object = lib2)
lib2 <- RunUMAP(object = lib2, dims = 1:50)




## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_lib2 <- paramSweep_v3(lib2, PCs = 1:50, sct = FALSE)
sweep.res.list_lib2
sweep.stats_lib2 <- summarizeSweep(sweep.res.list_lib2, GT = FALSE)
sweep.stats_lib2
bcmvn_lib2 <- find.pK(sweep.stats_lib2)

ggplot(bcmvn_lib2, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_lib2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))
pK  # 0.01

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- lib2@meta.data$seurat_clusters
annotations
homotypic.prop <- modelHomotypic(annotations)           
homotypic.prop
nExp_poi <- round(0.13*nrow(lib2@meta.data))  ## Assuming 13%% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
lib2 <- doubletFinder_v3(lib2,PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
lib2@meta.data
pdf("visualuzation_D_Slib2.pdf")
DimPlot(lib2, reduction = 'umap', group.by = "DF.classifications_0.25_0.01_1691")
dev.off()
pdf("dpn_distrubutionlib2.pdf")
FeaturePlot(lib2 , features = "LOC6627282")
dev.off()
# number of singlets and doublets
table(lib2@meta.data$ DF.classifications_0.25_0.01_1691)
#Doublet Singlet 
#1691   14536  

object_NoDoublets <- subset(lib2, subset =  DF.classifications_0.25_0.01_1691== "Singlet")
object_NoDoublets
saveRDS(object_NoDoublets , "lib2_NoDoublet.rds")
object_NoDoublets<- readRDS("lib2_NoDoublet.rds")

object_NoDoublets









