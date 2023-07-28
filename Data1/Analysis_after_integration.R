larvaOL.integrated <- readRDS("larvaOL.integrated_after_findNeighbourds_df.rds")
larvaOL.integrated
#larvaOL.integrated <-FindClusters(object = larvaOL.integrated, resolution = seq(0.5,3, by=0.5))
larvaOL.integrated <-FindClusters(object = larvaOL.integrated, resolution = 6 )

#pdf("clustree_integratedDatadf.pdf")
#clustree(larvaOL.integrated@meta.data[,grep("integrated_snn_res.", colnames(larvaOL.integrated@meta.data))], prefix = "integrated_snn_res.")
#dev.off()
#DimPlot(larvaOL.integrated, reduction = "umap",group.by = "integrated_snn_res.1.5", label = TRUE)


DefaultAssay(larvaOL.integrated) <-"RNA"
pdf(file = "Dimplotintegrateddf42.pdf")
DimPlot(larvaOL.integrated, reduction = "umap", label = TRUE)
dev.off()

pdf(file = "featureplotintegrateddf42.pdf")
FeaturePlot(larvaOL.integrated, features = "rna_LOC6627282" )
dev.off()

pdf(file = "vlnplotintegrateddf42.pdf", height = 8, width = 4)
VlnPlot(larvaOL.integrated, features = "LOC6627282", assay = "RNA", ncol = 2, pt.size = 0.1)
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
cluster_assignments <- larvaOL.integrated$integrated_snn_res.6[expressing_cells]
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
write.csv(cell_counts, file = "count_of_cell_intgration42.csv")
saveRDS(larvaOL.integrated , "larvaOL.integrated42.rds")
larvaOL.integrated42 <- readRDS("larvaOL.integrated42.rds")

###subset cluster of cells expressing dpn
nw_object <- c(42)
cluster42 <- subset(larvaOL.integrated42, idents = nw_object)
cluster42
saveRDS(cluster42, "cluster42.rds")

cluster42 <- readRDS("cluster42.rds")
FeaturePlot(cluster42 , features = c('LOC6629478' , 'LOC6629707') , blend = TRUE )
FeaturePlot(cluster42 , features = 'LOC6627282' )
FeaturePlot(cluster42 , features = c('LOC6629478' , 'LOC6629707'))
gene <- "LOC6627282"
expr <- FetchData(cluster42,vars = gene)
my_object <- cluster42[, which(x = expr > 0)]
my_object
FeaturePlot(my_object , features = c('LOC6629478' , 'LOC6629707') , blend = T )
FeaturePlot(my_object , features = 'LOC6627282' )
FeaturePlot(my_object , features = c('LOC6629478' , 'LOC6629707'))
saveRDS(my_object,"my_object_42.rds")
my_object <- readRDS("my_object_42.rds")
cluster42 <- readRDS("cluster42.rds")




library(slingshot)
library(pheatmap)
library(tradeSeq)


sce <- as.SingleCellExperiment(cluster42, assay = "RNA")

sce <- slingshot(sce, reducedDim = 'UMAP')
set.seed(3)
library(tradeSeq)
icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = sce@colData@listData[["slingshot"]]@assays@data@listData[["pseudotime"]],
                   cellWeights = sce@colData@listData[["slingshot"]]@assays@data@listData[["weights"]],
                   nGenes = 300,
                   k = 3:7)

set.seed(3)
sce <- fitGAM(sce,
              nknots = 7)
mean(rowData(sce)$tradeSeq$converged)

saveRDS(sce, "sce_cluster42.rds")

sce <- readRDS("sce_cluster42.rds")
sce
#run association test
ATres <- associationTest(sce)

saveRDS(ATres, "Atres_cluster42.rds")

ATres <- readRDS("Atres_cluster42.rds")
ATres

virilisTFs<-read.csv("virilis_list_TFs.csv")
virilisTFs_cand <- virilisTFs[,1]
is.vector(virilisTFs_cand)
head(virilisTFs_cand)

plotSmoothers(sce, assays(sce)$counts, gene = "LOC6636198", alpha =1, border = TRUE) + ggtitle("LOC6636198") #EY
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6629707", alpha =1, border = TRUE) + ggtitle("LOC6629707") #hth
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6629478", alpha =1, border = TRUE) + ggtitle("LOC6629478") #tll
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6627282", alpha =1, border = TRUE) + ggtitle("LOC6627282") #dpn
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6635714", alpha =1, border = TRUE) + ggtitle("LOC6635714") # slp1
dir.create("plots_TEMPORALdfDATA1")
dir.create("plots2_NONTEMPORALdfDATA1")

# Loop over the candidate TFs
for (i in 1:length(virilisTFs_cand)) {
  gene_id <- virilisTFs_cand[i]
  # Get the row metadata of the `sce` object
  rowData <- rowData(sce)
  
  # Extract the gene IDs from the row metadata
  gene_ids <- rownames(rowData)
  
  # Check whether the gene IDs in `virilisTFs_cand` are present in the row metadata
  if(gene_id %in% gene_ids){
    plot_data <- plotSmoothers(sce, assays(sce)$counts, gene = gene_id, alpha = 0.5, border = TRUE) + ggtitle(gene_id)
    x <- plot_data$data$time
    y <- plot_data$data$gene_count
    y_log <- log(plot_data$data$gene_count+1)
    x_max <- max(x)
    x_min <- min(x)
    # Divide the x axis into 3 intervals, try 6 intervals
    #num_intervals <- 4
    num_intervals <- 6
    interval_size <- (x_max - x_min) / num_intervals
    interval_starts <- seq(from = x_min, to = x_max -interval_size, by = interval_size)
    
    
    # Calculate the mean expression for each interval
    mean_vals <- tapply(y_log, findInterval(x, interval_starts), mean) # find intervalle containing each element of x in intervalle start
    # find the index position of x in xmin to intervallstart
    mean_vals
    mediane_vals <- tapply(y_log, findInterval(x, interval_starts), median)# apply function over subset of vector
    mediane_vals
    # vector/list of factor /fct
    
    # Find the index of the interval with the highest mean expression(location)
    max_index <- which.max(mean_vals)
    max_index
    
    min_index <- which.min(mean_vals)
    min_index
    
    if(mean_vals[max_index] > 0.5 && any(mediane_vals[-max_index] < 0.25) ) {
      ggsave(paste0("plots_TEMPORALdfDATA1/", gene_id, ".png"), plot = plot_data, width = 10, height = 7)
    }
    else{
      ggsave(paste0("plots2_NONTEMPORALdfDATA1/", gene_id, ".png"), plot = plot_data, width = 10, height = 7)
      
    }
    
  }
}

#####FetchData

nbs <- readRDS("my_object_42.rds")
nbs
sce <- as.SingleCellExperiment(nbs, assay = "RNA")

sce <- slingshot(sce, reducedDim = 'UMAP')
set.seed(3)
pdf("evaluateKnot")
icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = sce@colData@listData[["slingshot"]]@assays@data@listData[["pseudotime"]],
                   cellWeights = sce@colData@listData[["slingshot"]]@assays@data@listData[["weights"]],
                   nGenes = 300,
                   k = 3:7)
dev.off()
set.seed(3)
sce <- fitGAM(sce,
              nknots = 7)
mean(rowData(sce)$tradeSeq$converged)

saveRDS(sce, "sce_cluster42fetch.rds")

sce <- readRDS("sce_cluster42fetch.rds")
sce
#run association test
ATres <- associationTest(sce)

saveRDS(ATres, "Atres_cluster42fetch.rds")

ATres <- readRDS("Atres_cluster42fetch.rds")
ATres

virilisTFs<-read.csv("virilis_list_TFs.csv")
virilisTFs_cand <- virilisTFs[,1]
is.vector(virilisTFs_cand)
head(virilisTFs_cand)

plotSmoothers(sce, assays(sce)$counts, gene = "LOC6636198", alpha =1, border = TRUE) + ggtitle("LOC6636198") #EY
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6629707", alpha =1, border = TRUE) + ggtitle("LOC6629707") #hth
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6629478", alpha =1, border = TRUE) + ggtitle("LOC6629478") #tll
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6627282", alpha =1, border = TRUE) + ggtitle("LOC6627282") #dpn
plotSmoothers(sce, assays(sce)$counts, gene = "LOC6635714", alpha =1, border = TRUE) + ggtitle("LOC6635714") # slp1


dir.create("plots_TEMPORALdfDATA1fetch")
dir.create("plots2_NONTEMPORALdfDATA1fetch")

# Loop over the candidate TFs
for (i in 1:length(virilisTFs_cand)) {
  gene_id <- virilisTFs_cand[i]
  # Get the row metadata of the `sce` object
  rowData <- rowData(sce)
  
  # Extract the gene IDs from the row metadata
  gene_ids <- rownames(rowData)
  
  # Check whether the gene IDs in `virilisTFs_cand` are present in the row metadata
  if(gene_id %in% gene_ids){
    plot_data <- plotSmoothers(sce, assays(sce)$counts, gene = gene_id, alpha = 0.5, border = TRUE) + ggtitle(gene_id)
    x <- plot_data$data$time
    y <- plot_data$data$gene_count
    y_log <- log(plot_data$data$gene_count+1)
    x_max <- max(x)
    x_min <- min(x)
    # Divide the x axis into 3 intervals, try 6 intervals
    #num_intervals <- 4
    num_intervals <- 6
    interval_size <- (x_max - x_min) / num_intervals
    interval_starts <- seq(from = x_min, to = x_max -interval_size, by = interval_size)
    
    
    # Calculate the mean expression for each interval
    mean_vals <- tapply(y_log, findInterval(x, interval_starts), mean) # find intervalle containing each element of x in intervalle start
    # find the index position of x in xmin to intervallstart
    mean_vals
    mediane_vals <- tapply(y_log, findInterval(x, interval_starts), median)# apply function over subset of vector
    mediane_vals
    # vector/list of factor /fct
    
    # Find the index of the interval with the highest mean expression(location)
    max_index <- which.max(mean_vals)
    max_index
    
    min_index <- which.min(mean_vals)
    min_index
    
    if(mean_vals[max_index] > 0.5 && any(mediane_vals[-max_index] < 0.25) ) {
      ggsave(paste0("plots_TEMPORALdfDATA1fetch/", gene_id, ".png"), plot = plot_data, width = 10, height = 7)
    }
    else{
      ggsave(paste0("plots2_NONTEMPORALdfDATA1fetch/", gene_id, ".png"), plot = plot_data, width = 10, height = 7)
      
    }
    
  }
}


