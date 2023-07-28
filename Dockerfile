# Utilisez une image de base R
FROM r-base

# Installe les dépendances système nécessaires
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev

# Installe les packages R depuis CRAN
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'pheatmap', 'clustree' , 'dplyr' , 'patchwork'), repos = 'https://cran.rstudio.com/')"

# Installe le package Seurat depuis CRAN
RUN R -e "install.packages('Seurat', repos = 'https://cran.rstudio.com/')"

# Installe le package remotes depuis CRAN
RUN R -e "install.packages('remotes', repos = 'https://cran.rstudio.com/')"

# Installe le package SingleCellExperiment depuis Bioconductor
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install('SingleCellExperiment')"

# Installe les packages slingshot et tradeSeq depuis Bioconductor
RUN R -e "BiocManager::install(c('slingshot', 'tradeSeq'))"

# Installe le package DoubletFinder depuis GitHub
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"


## to rerun : 
#docker build -t nouhaylaerraysy/single_cell:v1.0 .
#docker push  nouhaylaerraysy/single_cell:v1.0 