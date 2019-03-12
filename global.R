#devtools packages

list.of.dev.packages <- c("plotly",
                          "crosstalk",
                          "DT")


#new.packages.dev <- list.of.dev.packages[!(list.of.dev.packages %in% installed.packages()[,"Package"])]

#if( "plotly" %in% new.packages.dev) devtools::install_github("ropensci/plotly", upgrade = "never")
#if( "crosstalk" %in% new.packages.dev) devtools::install_github("rstudio/crosstalk",force=TRUE, upgrade = "never")
#if( "DT" %in% new.packages.dev) devtools::install_github('rstudio/DT', upgrade = "never")


library(plotly)
library(crosstalk)
library(DT)


# CRAN packages
list.of.packages <- c("shiny",
                      "shinyBS",
                      "RcppArmadillo",
                      "GSEABase",
                      "shinyjs",
                      'RColorBrewer',
                      "stringr",
                      'formula.tools',
                      'data.table',
                      'fdrtool',
                      "VennDiagram",
                      "devtools",
                      'colorspace',
                      "officer",
                      "magrittr",
                      "openxlsx",
                      "ggrepel",
                      "V8",
                      "WGCNA",
                      'svglite',
                      "visNetwork",
                      "ggpubr",
                      "gplots",
                      "bindrcpp",
                      "pheatmap",
                      "purrr"
)

#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#if(length(new.packages)>0) install.packages(new.packages, dependencies = T)

# BioconductoR packages
list.of.bioc.packages<- c("rhdf5",
                          "DESeq2",
                          "IHW",
                          "tximport",
                          'clusterProfiler',
                          "org.Hs.eg.db",
                          "org.Mm.eg.db",
                          "org.Mmu.eg.db",
                          "sva",
                          "limma",
                          "geneplotter",
                          "rhdf5",
                          'biomaRt',
                          "AnnotationDbi", 
                          "impute", 
                          "GO.db", 
                          "preprocessCore",
                          "pcaGoPromoter",
                          "pcaGoPromoter.Mm.mm9",
                          "pcaGoPromoter.Hs.hg19",
                          "pathview",
                          "lpsymphony",
                          "officer"
)



#new.packages.bioc <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]

#source("https://bioconductor.org/biocLite.R")
#if(length(new.packages.bioc)>0) biocLite(new.packages.bioc,suppressUpdates=TRUE)

lapply(c(list.of.dev.packages,list.of.packages,list.of.bioc.packages), require, character.only = TRUE)



