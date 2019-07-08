#devtools packages

options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")

list.of.dev.packages <- c("plotly",
                          "crosstalk",
                          "DT")


new.packages.dev <- list.of.dev.packages[!(list.of.dev.packages %in% installed.packages()[,"Package"])]

if( "plotly" %in% new.packages.dev) devtools::install_github("ropensci/plotly", upgrade = "never")
if( "crosstalk" %in% new.packages.dev) devtools::install_github("rstudio/crosstalk",force=TRUE, upgrade = "never")
if( "DT" %in% new.packages.dev) devtools::install_github('rstudio/DT', upgrade = "never")
# 
# 
# library(plotly)
# library(crosstalk)
# library(DT)
# library(shinyBS)
# library(RcppArmadillo)
# 
# library(GSEABase)
# library(shinyjs)
# library(RColorBrewer)
# library(stringr)
# library(formula.tools)
# 
# library(data.table)
# library(crosstalk)
# library(VennDiagram)
# library(fdrtool)
# library(devtools)
# 
# library(colorspace)
# library(officer)
# library(magrittr)
# library(openxlsx)
# library(ggrepel)
# 
# library(V8)
# library(WGCNA)
# library(DT)
# library(svglite)
# library(visNetwork)
# 
# library(ggpubr)
# library(gplots)
# library(bindrcpp)
# library(pheatmap)
# library(purrr)

#CRAN packages
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

lapply(list.of.packages, require, character.only = TRUE)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages, dependencies = T)


library(rhdf5)
library(DESeq2)
library(IHW)
library(tximport)
library(clusterProfiler)

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Mmu.eg.db)
library(sva)
library(limma)

library(geneplotter)
library(biomaRt)
library(AnnotationDbi)
library(impute)
library(GO.db)

library(preprocessCore)
library(pcaGoPromoter)
library(pathview)
library(pcaGoPromoter.Mm.mm9)
library(pcaGoPromoter.Hs.hg19)
library(lpsymphony)
 
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
                          'biomaRt',
                          "AnnotationDbi",
                          "impute",
                          "GO.db",
                          "preprocessCore",
                          "pcaGoPromoter",
                          "pcaGoPromoter.Mm.mm9",
                          "pcaGoPromoter.Hs.hg19",
                          "pathview",
                          "lpsymphony"
)


new.packages.bioc <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]

install.packages("BiocManager")
if(length(new.packages.bioc)>0) BiocManager::install(new.packages.bioc,suppressUpdates=TRUE)

lapply(c(list.of.dev.packages,list.of.packages,list.of.bioc.packages), require, character.only = TRUE)



