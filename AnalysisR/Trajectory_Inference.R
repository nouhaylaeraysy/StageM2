suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment)
  library(RColorBrewer); library(scales)
  library(viridis)
  library(pheatmap)
  library(knitr); library(gridExtra)
  library(tradeSeq); library(RColorBrewer)
})

cluster_subset <- readRDS("clusterNBS.rds")
cluster_subset
sce <- as.SingleCellExperiment(cluster_subset, assay = "RNA")

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
              nknots = 6)
mean(rowData(sce)$tradeSeq$converged)

saveRDS(sce, "sce_clusterNBs.rds")

sce <- readRDS("sce_clusterNBs.rds")
sce
#run association test
ATres <- associationTest(sce)

saveRDS(ATres, "Atres_clusterNBs.rds")

ATres <- readRDS("Atres_clusterNBs.rds")
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


