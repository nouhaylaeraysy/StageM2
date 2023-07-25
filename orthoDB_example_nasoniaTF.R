library(httr)
library(jsonlite)
library(readr)

# list of genes (example)
gene_ids <- c('LOC100678172' ) #, 'LOC103315530' , 'LOC100678407')


# Loop to execute the search for each gene ID
for (gene_id in gene_ids) {
  url <- paste0("https://data.orthodb.org/v11/genesearch?query=", gene_id)
  
  # Send GET request and retrieve data
  response <- httr::GET(url)
  data <- httr::content(response, "text")
  
  # Store data in a file  
  write(data, file = paste0("result_", trimws(gene_id), ".txt"))  # remove space 
  
  cat(paste("Data for gene ID", gene_id, "has been saved in the file 'result_", gene_id, ".txt'.\n"))
  
  # Read the file
  data <- jsonlite::fromJSON(txt = readLines(paste0("result_", trimws(gene_id), ".txt")), simplifyVector = FALSE)
  
  # Extract orthologue  Drosophila melanogaster
  if ("orthologs_in_model_organisms" %in% names(data)) {
    orthologs <- data$orthologs_in_model_organisms
    if (!is.null(orthologs) && length(orthologs) > 0) {
      # Filter orthologs to get the one with Drosophila melanogaster
      drosophila_ortholog <- lapply(orthologs, function(ortholog) {  #filter out the orthologs that match the organism name "Drosophila melanogaster". 
        if (grepl("Drosophila melanogaster", ortholog$organism$name, ignore.case = TRUE)) {
          return(ortholog)
        }
      })
      
      drosophila_ortholog <- drosophila_ortholog[!sapply(drosophila_ortholog, is.null)] # removing any element that has value null
      
      if (length(drosophila_ortholog) > 0) {
        gene_id <- drosophila_ortholog[[1]]$genes[[1]]$gene_id$id
        extracted_value <- strsplit(gene_id, ";")[[1]][1]
        cat("Extracted:", extracted_value, "\n")
      } else {
        cat("No ortholog found for Drosophila melanogaster.\n")
      }
    } else {
      cat("No orthologs found for the gene ID.\n")
    }
  }
  
  # Wait for 5 seconds before the next request
  Sys.sleep(3)
}

