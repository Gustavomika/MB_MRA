# Medulloblastoma Master Regulator Analysis using RTN R/Bioconductor package.

This repository contains the scripts used to infer the medulloblastoma (MB) regulatory network and identification of master regulator for MB (subgroups SHH, group 3, group 4) and risk master regulators. Briefly, the network was infered using GSE85217 MB samples. Gene signitures were identified comparing a subset of GSE85217 tumor samples against GSE167447 healthy fetal cerebellum samples. Than, the MB network and gene signatures were used to identify the master regulator. Clinical data were used o identify we used risk master regulators.

Data used for this analysis are publicly available in Gene Expression Omnibus (GEO-NCBI). Before start, download than and extract at:

"./Dados/Brutos/Amostras/GSE167447" for GSE167447 samples;
"./Dados/Brutos/Amostras/GSE85217" for GSE85217 samples.

Links to data specified on the scripts: 
* [GSE85217](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85217)
* [GSE167447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167447)

For package installing, there is two sources: The Comprehensive R Archive Network - [CRAN](https://cran.r-project.org/) - and [Bioconductor](https://bioconductor.org). To install CRAN packages, on the R console, just type:

```{r}
install.packages("package")
```

For Bioconductor packages, the first installation requires the package BiocManager (from CRAN). After installing `BiocManager`, you can install Bioconductor packages using the following command on R console:

```{r}
BiocManager::install("package")
```
The parameter `package` can be substituted by an character vector (`c()`) containing all packages to install. 

R Packages used in this analysis:
1. From CRAN
* BiocManager
* data.table
* snow
* RColorBrewer 
* broom
* ggplot2 
* dplyr
* maxstat
* RColorBrewer
* survminer
* ggpubr
* openxlsx
* ggstatsplot
* survival
* pheatmap

2. From Bioconductor
* oligo
* limma
* arrayQualityMetrics
* AnnotationDbi
* RTN
* RedeR 
* clusterProfiler
* org.Hs.eg.db
* ComplexHeatmap
* RTNsurvival

