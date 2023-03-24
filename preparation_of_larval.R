library(Seurat)
library(patchwork)
library(clustree)
library(tradeSeq)

# Load the Larval data 
larval1.data <- Read10X(data.dir = "/home/adminkonstantinides9/Documents/for_nouhayla/larval1/filtered_feature_bc_matrix/")
larval2.data <- Read10X(data.dir = "/home/adminkonstantinides9/Documents/for_nouhayla/larval2/filtered_feature_bc_matrix/")
larval3.data <- Read10X(data.dir = "/home/adminkonstantinides9/Documents/for_nouhayla/larval3/filtered_feature_bc_matrix/")
dim(larval1.data)
#[1] 16102 13269  gene * cell
# Create seurat object with some low filtering: genes expressed in at least 3 cells, cells with at least 200 genes:
larval1 <- CreateSeuratObject(counts = larval1.data, project = "larval1", min.cells = 3 , min.features = 200)
larval1
larval2 <- CreateSeuratObject(counts = larval2.data, project = "larval1", min.cells = 3 , min.features = 200)
larval2
larval3 <- CreateSeuratObject(counts = larval3.data, project = "larval1", min.cells = 3 , min.features = 200)
larval3
# standard pre-processing workflow 

# step 1 : Quality controle and selecting cell for further analysis : 
# Visualize QC metrics as a violin plot : plot of number of counts and detected genes per sample
VlnPlot(larval1, features = c("nFeature_RNA", "nCount_RNA"))
VlnPlot(larval2, features = c("nFeature_RNA", "nCount_RNA"))
VlnPlot(larval3, features = c("nFeature_RNA", "nCount_RNA"))

# thresholds to filter out cells , based on vlnplot
larval1 <- subset(larval1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
larval1
#An object of class Seurat 
#12083 features across 13259 samples within 1 assay 
#Active assay: RNA (12083 features, 0 variable features)

larval2 <- subset(larval2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
larval2
larval3 <- subset(larval3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
larval3
#step2 : Normalization and scaling:
larval1 <- NormalizeData(larval1,
              normalization.method = "LogNormalize",
              scale.factor = 10000) # logarithme of gene expression for each cell *scaling factor /total expression 
larval2 <- NormalizeData(larval2,
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
larval3 <- NormalizeData(larval3,
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
#chek data after normalization 
GetAssayData(larval1)[1:10,1:10]
GetAssayData(larval2)[1:10,1:10]
GetAssayData(larval3)[1:10,1:10]

#step3 : Identification of highly variable features
#gene that are highly expressed in some cells and less in other 
larval1 <- FindVariableFeatures(larval1,selection.method = "vst", nfeatures = 2000)
larval2 <- FindVariableFeatures(larval1,selection.method = "vst", nfeatures = 2000)
larval3 <- FindVariableFeatures(larval1,selection.method = "vst", nfeatures = 2000)
# top 10 most variable genes
top10 <- head(VariableFeatures(larval1), 10)
top10
# [1] "LOC6630360"   "LOC6632158"   "LOC116651881" "LOC26530589"  "LOC6633783"   "LOC6635453"   "LOC6635718"   "LOC6622406"   "LOC6628401"   "LOC6632157"  
# we have the same genes in the larval2 and larval3 
# plot as scatterplot of average expression level (x-axis) versus variance (y-axis)
# store unlabeled plot in object vf_plot
vf_plot <- VariableFeaturePlot(larval1) 
vf_plot
# add gene labels of top variable genes
LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE) 

#step4:  scaling for PCA, only the variable genes are scaled by default : linear transformation with the conservation of biological differences
larval1 <- ScaleData(larval1, verbose = FALSE) 
larval2 <- ScaleData(larval2, verbose = FALSE) 
larval3 <- ScaleData(larval3, verbose = FALSE) 

#step5 : Perform linear dimensional reduction
#Run PCA and plot 
larval1 <- RunPCA(larval1, npcs = 200, verbose = FALSE) 
DimPlot(larval1, reduction = "pca")
larval2<- RunPCA(larval2, npcs = 200, verbose = FALSE) 
DimPlot(larval2, reduction = "pca")
larval3<- RunPCA(larval3, npcs = 200, verbose = FALSE) 
DimPlot(larval3, reduction = "pca")

# Heatmap of gene expression of top genes contributing to each of the 12 first PCs:
DimHeatmap(larval1, dims = 1:50, cells = 500, balanced = TRUE)
DimHeatmap(larval1, dims = 1:20, cells = 250, balanced = TRUE)
DimHeatmap(larval1, dims = 21:40, cells = 500, balanced = TRUE)
DimHeatmap(larval1, dims = 41:50, cells = 500, balanced = TRUE)
DimHeatmap(larval1, dims = 51:60, cells = 250, balanced = TRUE)
#Elbowplot
ElbowPlot(larval1, ndims = 150)
ElbowPlot(larval1, ndims = 200)
 
 # step6 : cluster the cells : 
# 1 -larval1
larval1 <- FindNeighbors(larval1, dims = 1:135)
larval1 <- FindClusters(larval1, resolution = seq(0.5,3, by=0.5))
clustree(larval1@meta.data[,grep("RNA_snn_res.", colnames(larval1@meta.data))], prefix = "RNA_snn_res.")
DimPlot(larval1, group.by = "RNA_snn_res.1.5")
FeaturePlot(larval1, features = "LOC6627282" ) # check deadpan expression
larval1 <- RunUMAP(larval1,dims = 1:135) # with 135 dimension
DimPlot(larval1, reduction = "umap")
FeaturePlot(larval1, features = "LOC6627282" )
VlnPlot(larval1, features = "LOC6627282",assay = "RNA", ncol = 2, pt.size = 0.1, group.by = "RNA_snn_res.1.5")
### e.g DEADPAN
deadpan <- "LOC6627282"

# Subset the data to only include cells expressing the gene of interest
gene_expression <- larval1@assays$RNA@data[deadpan, ]
head(gene_expression)
# Get the indices of cells that express the gene of interest
expressing_cells <- which(gene_expression > 0)
head(expressing_cells)
# Subset the cluster assignments of the expressing cells
cluster_assignments <- larval1$RNA_snn_res.2[expressing_cells]
head(cluster_assignments)
# Subset the original identifications of the expressing cells
cell_idents <- larval1$orig.ident[expressing_cells]
cell_idents
# Combine the cluster assignments and cell idents into a data frame
cell_data <- data.frame(Cluster = cluster_assignments, ID = cell_idents)
head(cell_data)
# Get the counts of cells per cluster
cell_counts <- table(cell_data)

# Print the table of cluster sizes
print(cell_counts)

#2 - larval2
larval2 <- FindNeighbors(larval2, dims = 1:135)
larval2 <- FindClusters(larval2, resolution = seq(0.5,3, by=0.5))
clustree(larval2@meta.data[,grep("RNA_snn_res.", colnames(larval2@meta.data))], prefix = "RNA_snn_res.")
DimPlot(larval2, group.by = "RNA_snn_res.1.5")
larval2 <- RunUMAP(larval2,dims = 1:135) # with 135 dimension
DimPlot(larval2, reduction = "umap")
FeaturePlot(larval2, features = "LOC6627282" )
VlnPlot(larval2, features = "LOC6627282",assay = "RNA", ncol = 2, pt.size = 0.1, group.by = "RNA_snn_res.1.5")

gene_expression <- larval2@assays$RNA@data[deadpan, ]
expressing_cells <- which(gene_expression > 0)

# Subset the cluster assignments of the expressing cells
cluster_assignments <- larval2$RNA_snn_res.2[expressing_cells]
head(cluster_assignments)

# Subset the original identifications of the expressing cells
cell_idents <- larval2$orig.ident[expressing_cells]

# Combine the cluster assignments and cell idents into a data frame
cell_data <- data.frame(Cluster = cluster_assignments, ID = cell_idents)

# Get the counts of cells per cluster
cell_counts <- table(cell_data)

# Print the table of cluster sizes
print(cell_counts)


# 3- larval3
larval3 <- FindNeighbors(larval3, dims = 1:135)
larval3 <- FindClusters(larval3, resolution = seq(0.5,3, by=0.5))
clustree(larval3@meta.data[,grep("RNA_snn_res.", colnames(larval3@meta.data))], prefix = "RNA_snn_res.")
DimPlot(larval3, group.by = "RNA_snn_res.1.5")
FeaturePlot(larval3, features = "LOC6627282" )
larval3 <- RunUMAP(larval1,dims = 1:135) # with 135 dimension
DimPlot(larval3, reduction = "umap")
FeaturePlot(larval3, features = "LOC6627282" )
VlnPlot(larval3, features = "LOC6627282",assay = "RNA", ncol = 2, pt.size = 0.1, group.by = "RNA_snn_res.1.5")


gene_expression <- larval3@assays$RNA@data[deadpan, ]
expressing_cells <- which(gene_expression > 0)

# Subset the cluster assignments of the expressing cells
cluster_assignments <- larval3$RNA_snn_res.1.5[expressing_cells]
head(cluster_assignments)

# Subset the original identifications of the expressing cells
cell_idents <- larval3$orig.ident[expressing_cells]

# Combine the cluster assignments and cell idents into a data frame
cell_data <- data.frame(Cluster = cluster_assignments, ID = cell_idents)

# Get the counts of cells per cluster
cell_counts <- table(cell_data)

# Print the table of cluster sizes
print(cell_counts)
