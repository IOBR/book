
# **How to install IOBR**

## ⏳ Install Dependency Packages

It is essential that you have R 3.6.3 or above already installed on your computer or server. IOBR is a pipeline that utilizes many other R packages that are currently available from CRAN, Bioconductor and GitHub.


```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
          "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
          "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr')
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}
```

## 💻 Install IOBR package

When the dependent environments are built, users are able to install IOBR from github by typing the following code into your R session:

```r
if (!requireNamespace("IOBR", quietly = TRUE))  devtools::install_github("IOBR/IOBR")

library(IOBR)
```


## 📌 How to update IOBR

```r
detach("package:IOBR")
path<-.libPaths()
remove.packages(c('IOBR'), lib=file.path(path))
devtools::install_github("IOBR/IOBR")
```

## 🔍 The main pipeline of IOBR

<div class="figure" style="text-align: center">
<img src="./fig/IOBR-Package.png" alt="The main pipeline of IOBR" width="95%" />
<p class="caption">(\#fig:flowchart)The main pipeline of IOBR</p>
</div>

## ☝️ Main Functions of IOBR

* <div style="color:green">**Data Preparation: data annotation and transformation**</div> 
   * `count2tpm()`: transform count data of RNA sequencing into TPM data.
   * `anno_eset()`: annotate the normalized genes expression matrix, including RNAseq and array (Affymetrix or Illumina).
   * `remove_duplicate_genes()`: remove the genes annotated with the duplicated symbol after normalization and retain only the symbol with highest expression level.
   * `mouse2human_eset()`: Converting muouse gene symbol to human gene symbol of expression set.
   * `find_outlier_samples()`: Waiting for updates...
   * `remove_batcheffect()`: Waiting for updates...
   </br>

* <div style="color:green">**TME Deconvolution Module: integrate multiple algorithms to decode immune contexture**</div> 
   * `deconvo_tme()`: decode the TME infiltration with different deconvolution methodologies, based on bulk RNAseq, microarray or single cell RNAseq data.
   * `generateRef()`: generate a novel gene reference matrix for a specific feature such as infiltrating cell, through the  SVR and lsei algorithm.
</br>

* <div style="color:green">**Signature Module: calculatint signature scores, estimate phenotype related signatures and corresponding genes, and evaluate signatures generated from single-cell RNA sequencing data **</div>
  * `calculate_sig_score()`: estimate the interested signatures enrolled in IOBR R package, which involves TME-associated, tumor-metabolism, and tumor-intrinsic signatures.
  * `feature_manipulation()`: manipulate features including the cell fraction and signatures generated from multi-omics data for latter analysis and model construction. Remove missing values, outliers and variables without significant variance.
  * `format_signatures()`: generate the object for `calculate_sig_score()` function, by inputting a data frame with signatures as column names of corresponding gene sets, and return a list contain the signature information for calculating multiple signature scores.
  * `format_msigdb()`: transform the signature gene sets data  with gmt format, which is not included in the signature collection and might be downloaded in the MSgiDB website, into the object of `calculate_sig_score()`function.
  * `sig_gsea()`: Waiting for updates...
</br>
    
* <div style="color:green">**Batch Analysis and Visualization: batch survival analysis and batch correlation analysis and other batch statistical analyses **</div>
  * `batch_surv()`: batch survival analysis of multiple continuous variables including varied signature scores.
  * `subgroup_survival()`: batch survival analysis of multiple categorized variables with different number of subgroups. 
  * `batch_cor()`: batch analysis of correlation between two continuous variables using Pearson correlation coefficient or Spearman's rank correlation coefficient .
  * `batch_wilcoxon()`: conduct batch wilcoxon analyses of binary variables.
  * `batch_pcc()`: batch analyses of Partial Correlation coefficient(PCC) between continuous variables  and minimize the interference derived from confounding factors.
  * `iobr_cor_plot()`: visualization of batch correlation analysis of signatures from 'sig_group'. Visualize the correlation between signature or phenotype with  expression of gene sets  in target signature is also supported.
  * `cell_bar_plot()`: batch visualization of TME cell fraction, supporting input of deconvolution results from 'CIBERSORT', 'EPIC' and 'quanTIseq' methodologies to further compare the TME cell distributions within one sample or among different samples.
  * `iobr_pca()`: The iobr_pca function performs Principal Component Analysis (PCA), which reduces the dimensionality of data while maintaining most of the original variance, and visualizes the PCA results on a scatter plot.
  * `iobr_cor_plot()`: Integrative correlation between phenotype and features.
  * `iobr_deg()`: Waiting for updates...
  * `get_cor()`: Waiting for updates...
  * `roc_time()`: Waiting for updates...
  * `sig_box()`: Waiting for updates...
  * `sig_heatmap()`: Waiting for updates...
  * `sig_forest()`: Waiting for updates...
  * `sig_roc()`: Waiting for updates...
  * `sig_surv_plot()`: Waiting for updates...
  * `find_markers_in_bulk()`: Waiting for updates...
  
</br>

* <div style="color:green">**Signature Associated Mutation Module: identify and analyze mutations relevant to targeted signatures**</div>
  * `make_mut_matrix()`: transform the mutation data with MAF format(contain the columns of gene ID and the corresponding gene alterations which including SNP, indel and frameshift) into a mutation matrix in a suitable manner for further investigating signature relevant mutations.
  * `find_mutations()`: identify mutations associated with a distinct phenotype or signature.
</br>

* <div style="color:green">**Model Construction Module: feature selection and fast model construct to predict clinical phenotype**</div>
  * `BinomialModel()`: select features and construct a model to predict a binary phenotype.
  * `PrognosticMode()`: select features and construct a model to predict clinical survival outcome.
</br>


## 🌐 Current working environment


```r
sessionInfo()
```

```
## R version 4.2.0 (2022-04-22 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19045)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
## [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
## [3] LC_MONETARY=Chinese (Simplified)_China.utf8
## [4] LC_NUMERIC=C                               
## [5] LC_TIME=Chinese (Simplified)_China.utf8    
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] IOBR_0.99.9           survminer_0.4.9       patchwork_1.1.1      
##  [4] clusterProfiler_4.4.3 survival_3.3-1        ggpubr_0.4.0         
##  [7] ggplot2_3.4.2         dplyr_1.1.2           tibble_3.2.1         
## [10] tidyHeatmap_1.8.1     ComplexHeatmap_2.15.4
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                  tidyselect_1.2.0           
##   [3] RSQLite_2.2.14              AnnotationDbi_1.58.0       
##   [5] BiocParallel_1.30.3         lpSolve_5.6.15             
##   [7] scatterpie_0.1.7            ScaledMatrix_1.4.0         
##   [9] munsell_0.5.0               codetools_0.2-18           
##  [11] preprocessCore_1.58.0       withr_2.5.0                
##  [13] colorspace_2.0-3            GOSemSim_2.22.0            
##  [15] Biobase_2.56.0              limSolve_1.5.6             
##  [17] highr_0.9                   knitr_1.39                 
##  [19] rstudioapi_0.13             SingleCellExperiment_1.18.0
##  [21] stats4_4.2.0                ggsignif_0.6.3             
##  [23] DOSE_3.22.0                 MatrixGenerics_1.8.0       
##  [25] GenomeInfoDbData_1.2.8      KMsurv_0.1-5               
##  [27] polyclip_1.10-0             bit64_4.0.5                
##  [29] farver_2.1.0                rhdf5_2.40.0               
##  [31] downloader_0.4              treeio_1.21.2              
##  [33] vctrs_0.6.2                 generics_0.1.3             
##  [35] xfun_0.40                   R6_2.5.1                   
##  [37] doParallel_1.0.17           GenomeInfoDb_1.34.4        
##  [39] clue_0.3-61                 graphlayouts_0.8.0         
##  [41] rsvd_1.0.5                  locfit_1.5-9.5             
##  [43] rhdf5filters_1.8.0          gridGraphics_0.5-1         
##  [45] bitops_1.0-7                cachem_1.0.6               
##  [47] fgsea_1.22.0                DelayedArray_0.22.0        
##  [49] assertthat_0.2.1            scales_1.2.0               
##  [51] ggraph_2.0.5                enrichplot_1.16.1          
##  [53] googlesheets4_1.0.0         gtable_0.3.1               
##  [55] beachmat_2.12.0             tidygraph_1.2.1            
##  [57] rlang_1.1.0                 genefilter_1.78.0          
##  [59] GlobalOptions_0.1.2         splines_4.2.0              
##  [61] lazyeval_0.2.2              rstatix_0.7.0              
##  [63] gargle_1.2.0                broom_0.8.0                
##  [65] yaml_2.3.5                  reshape2_1.4.4             
##  [67] abind_1.4-5                 modelr_0.1.8               
##  [69] backports_1.4.1             qvalue_2.28.0              
##  [71] tools_4.2.0                 bookdown_0.35              
##  [73] ggplotify_0.1.0             ellipsis_0.3.2             
##  [75] RColorBrewer_1.1-3          proxy_0.4-27               
##  [77] BiocGenerics_0.42.0         Rcpp_1.0.9                 
##  [79] plyr_1.8.7                  sparseMatrixStats_1.8.0    
##  [81] zlibbioc_1.42.0             purrr_0.3.4                
##  [83] RCurl_1.98-1.7              GetoptLong_1.0.5           
##  [85] viridis_0.6.2               cowplot_1.1.1              
##  [87] S4Vectors_0.34.0            zoo_1.8-10                 
##  [89] SummarizedExperiment_1.26.1 haven_2.5.0                
##  [91] ggrepel_0.9.1               cluster_2.1.3              
##  [93] fs_1.5.2                    magrittr_2.0.3             
##  [95] data.table_1.14.2           DO.db_2.9                  
##  [97] circlize_0.4.15             reprex_2.0.1               
##  [99] googledrive_2.0.0           matrixStats_0.62.0         
## [101] GSVA_1.44.2                 hms_1.1.1                  
## [103] evaluate_0.15               xtable_1.8-4               
## [105] XML_3.99-0.10               readxl_1.4.0               
## [107] IRanges_2.30.0              gridExtra_2.3              
## [109] shape_1.4.6                 compiler_4.2.0             
## [111] shadowtext_0.1.2            crayon_1.5.2               
## [113] htmltools_0.5.2             ggfun_0.0.6                
## [115] tzdb_0.3.0                  tidyr_1.2.0                
## [117] geneplotter_1.74.0          aplot_0.1.6                
## [119] lubridate_1.8.0             DBI_1.1.2                  
## [121] tweenr_1.0.2                corrplot_0.92              
## [123] dbplyr_2.2.0                MASS_7.3-56                
## [125] Matrix_1.5-4.1              car_3.1-0                  
## [127] readr_2.1.2                 cli_3.4.1                  
## [129] quadprog_1.5-8              parallel_4.2.0             
## [131] igraph_1.3.2                GenomicRanges_1.48.0       
## [133] forcats_0.5.1               pkgconfig_2.0.3            
## [135] km.ci_0.5-6                 xml2_1.3.3                 
## [137] foreach_1.5.2               ggtree_3.4.4               
## [139] annotate_1.74.0             XVector_0.36.0             
## [141] rvest_1.0.2                 yulab.utils_0.0.4          
## [143] stringr_1.4.0               digest_0.6.29              
## [145] graph_1.74.0                Biostrings_2.64.0          
## [147] rmarkdown_2.14              cellranger_1.1.0           
## [149] fastmatch_1.1-3             tidytree_0.3.9             
## [151] survMisc_0.5.6              dendextend_1.15.2          
## [153] DelayedMatrixStats_1.18.0   GSEABase_1.58.0            
## [155] rjson_0.2.21                nlme_3.1-157               
## [157] lifecycle_1.0.3             jsonlite_1.8.0             
## [159] Rhdf5lib_1.18.2             carData_3.0-5              
## [161] viridisLite_0.4.1           limma_3.52.1               
## [163] fansi_1.0.3                 pillar_1.9.0               
## [165] lattice_0.20-45             KEGGREST_1.36.2            
## [167] fastmap_1.1.0               httr_1.4.3                 
## [169] GO.db_3.15.0                emo_0.0.0.9000             
## [171] glue_1.6.2                  png_0.1-7                  
## [173] iterators_1.0.14            glmnet_4.1-4               
## [175] bit_4.0.4                   HDF5Array_1.24.1           
## [177] ggforce_0.3.3               class_7.3-20               
## [179] stringi_1.7.6               blob_1.2.3                 
## [181] BiocSingular_1.12.0         DESeq2_1.36.0              
## [183] memoise_2.0.1               irlba_2.3.5                
## [185] ape_5.6-2                   tidyverse_1.3.2            
## [187] e1071_1.7-11
```
