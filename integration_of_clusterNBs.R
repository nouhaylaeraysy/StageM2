library(seurat)
NBS1 <- -readRDS("~/clusterNBs1.rds")
NBS2 <- -readRDS("~/clusterNBs2.rds")
NBS3 <- -readRDS("~/clusterNBs3.rds")

reference_list <- c(NBS1,NBS2, NBS3)
larvalOl.anchors <- FindIntegrationAnchors(object.list = reference_list , dims = 1:150)
larvalOl.integrated <- IntegrateData(anchorset = larvalOl.anchors)
larvalOl.integrated <- FindVariableFeatures(larvalOl.integrated, assay = "RNA",
                            selection.method = "vst",
                            nfeatures = 2000,
                            verbose = FALSE)
larvalOl.integrated <- ScaleData(object = larvalOl.integrated , verbose = FALSE)
larvalOl.integrated <- RunPCA(object = larvalOl.integrated , npcs = 200 , verbose = FALSE)
larvalOl.integrated <- RunUMAP(object = larvalOl.integrated , dims = 1:150)
DimPlot(larvalOl.integrated, reduction = "umap")
larvalOl.integrated <- FindNeighbors(object = larvalOl.integrated , dims = 1:150)
larvalOl.integrated <- FindClusters(object = larvalOl.integrated , resolution = 2)
DimPlot(larvalOl.integrated, label = TRUE) + NoLegend()


FeaturePlot(larvalOl.integrated, features = "LOC6627282")
