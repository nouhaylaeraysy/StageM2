#Load libraries
library(Seurat) # 4.3.0
library(patchwork)
library(clustree)

# Load the data 
lib1 <- readRDS("lib1_NoDoublet.rds")
lib2 <- readRDS("lib2_NoDoublet.rds")
lib3 <- readRDS("lib3_NoDoublet.rds")
reference_list <- c(lib1,lib2, lib3)
reference_list

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = reference_list)
larvaOL.anchors <- FindIntegrationAnchors(object.list = reference_list, dims = 1:150 , anchor.features = features)
saveRDS(larvaOL.anchors ,"larvaOL.anchors_doubletFinder.rds")
larvaOL.anchors <- readRDS("larvaOL.anchors_doubletFinder.rds")
larvaOL.integrated <- IntegrateData(anchorset = larvaOL.anchors, dims = 1:150)
saveRDS(larvaOL.integrated , "larvalolintegrated_from_optimlib_doubletFinder.rds")
larvaOL.integrated <- readRDS("larvalolintegrated_from_optimlib_doubletFinder.rds")
DefaultAssay(larvaOL.integrated) <- "integrated"
larvaOL.integrated <- ScaleData(object = larvaOL.integrated, verbose = FALSE)
larvaOL.integrated <- RunPCA(object = larvaOL.integrated, npcs = 200, verbose = FALSE)
ElbowPlot(larvaOL.integrated, ndims = 200)
larvaOL.integrated <- RunUMAP(object = larvaOL.integrated, dims = 1:150)
larvaOL.integrated
pdf("dimplot_orig.identsDF.pdf")
DimPlot(larvaOL.integrated, reduction = "umap" , group.by = "orig.ident") # cells from both library are over each other are not separated due to conditions exactly what we want to doo with integration 
dev.off()
larvaOL.integrated
larvaOL.integrated <-FindNeighbors(object = larvaOL.integrated, dims = 1:150)
saveRDS(larvaOL.integrated , "larvaOL.integrated_after_findNeighbourds_df.rds")
larvaOL.integrated <- readRDS("larvaOL.integrated_after_findNeighbourds_df.rds")
#larvaOL.integrated <-FindClusters(object = larvaOL.integrated, resolution = seq(0.5,10, by=0.5))
larvaOL.integrated <-FindClusters(object = larvaOL.integrated, resolution = 1.5)

#pdf("clustree_integratedDatadf.pdf")
#clustree(larvaOL.integrated@meta.data[,grep("integrated_snn_res.", colnames(larvaOL.integrated@meta.data))], prefix = "integrated_snn_res.")
#dev.off()


DefaultAssay(larvaOL.integrated) <-"RNA"
pdf(file = "Dimplotintegrateddf.pdf")
DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE)
dev.off()

pdf(file = "featureplotintegrateddf.pdf")
FeaturePlot(larvaOL.integrated, features = "rna_LOC6627282" )
dev.off()

pdf(file = "vlnplotintegrateddf.pdf")
VlnPlot(larvaOL.integrated, features = "LOC6627282",assay = "RNA", ncol = 2, pt.size = 0.1)
dev.off()
## CLUSTER 2


### e.g DEADPAN
deadpan <- "LOC6627282"

# Subset the data to only include cells expressing the gene of interest
gene_expression <- larvaOL.integrated@assays$RNA@data[deadpan, ]
head(gene_expression)
# Get the indices of cells that express the gene of interest
expressing_cells <- which(gene_expression > 0)
head(expressing_cells)
# Subset the cluster assignments of the expressing cells
cluster_assignments <- larvaOL.integrated$integrated_snn_res.1.5[expressing_cells]
head(cluster_assignments)
# Subset the original identifications of the expressing cells
cell_idents <- larvaOL.integrated$orig.ident[expressing_cells]
cell_idents
# Combine the cluster assignments and cell idents into a data frame
cell_data <- data.frame(Cluster = cluster_assignments, ID = cell_idents)
head(cell_data)
# Get the counts of cells per cluster
cell_counts <- table(cell_data)

# Print the table of cluster sizes
print(cell_counts)  # check the cluster finded by VlnPlot
write.csv(cell_counts, file = "count_of_cell_intgration.csv")

###subset cluster of cells expressing dpn
nw_object <- c(29)
cluster29NBS <- subset(larvaOL.integrated, idents = nw_object)
cluster29NBS
saveRDS(cluster29NBS, "cluster29NBS.rds")

FeaturePlot(cluster29NBS , features = c('LOC6629478' , 'LOC6629707') , blend = TRUE )
FeaturePlot(cluster29NBS , features = 'rna_LOC6627282' )

gene <- "LOC6627282"
expr <- FetchData(cluster29NBS,vars = gene)
my_object <- cluster29NBS[, which(x = expr > 0)]
my_object
FeaturePlot(my_object , features = c('LOC6629478' , 'LOC6629707') , blend = T )
FeaturePlot(my_object , features = 'LOC6627282' )
FeaturePlot(my_object , features = c('LOC6629478' , 'LOC6629707'))
saveRDS(my_object,"my_object_FetchDataDa1.rds")
my_object <- readRDS("my_object_FetchDataDa1.rds")


