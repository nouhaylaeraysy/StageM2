library(httr)
library(jsonlite)
library(readr)

# Read gene IDs from the text file
gene_ids <- read_lines("list_TF_PheatmapKons1.txt")

# Loop to execute the search for each gene ID
for (gene_id in gene_ids) {
  url <- paste0("https://data.orthodb.org/v11/genesearch?query=", gene_id)
  
  # Send GET request and retrieve data
  response <- GET(url)
  data <- content(response, "text")
  
  # Store data in a file  
  write(data, file = paste0("result_", trimws(gene_id), ".txt"))
  
  cat(paste("Data for gene ID", gene_id, "has been saved in the file 'result_", gene_id, ".txt'.\n"))
  
  # Read the file
  # Read the file
  data <- fromJSON(txt = readLines(paste0("result_", trimws(gene_id), ".txt")), simplifyVector = FALSE)
  
  # Extract "CG8159" from "CG8159;Dmel\CG8159"
  if ("orthologs_in_model_organisms" %in% names(data)) {
    orthologs <- data$orthologs_in_model_organisms
    if (!is.null(orthologs) && length(orthologs) > 0) {
      gene_id <- orthologs[[1]]$genes[[1]]$gene_id$id
      extracted_value <- strsplit(gene_id, ";")[[1]][1]
      cat("Extracted:", extracted_value, "\n")
    }
  }
  
  # Wait for 5 seconds before the next request
  Sys.sleep(3)
}
# 