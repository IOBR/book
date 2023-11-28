
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

## 📇 Main Functions of IOBR

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


## 🌎 Current working environment


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
##  [1] IOBR_0.99.9           survminer_0.4.9       patchwork_1.1.3      
##  [4] survival_3.3-1        ggpubr_0.6.0          ggplot2_3.4.4        
##  [7] dplyr_1.1.2           tibble_3.2.1          clusterProfiler_4.4.4
## [10] tidyHeatmap_1.8.1     ComplexHeatmap_2.15.4
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                  tidyselect_1.2.0           
##   [3] RSQLite_2.3.2               AnnotationDbi_1.58.0       
##   [5] BiocParallel_1.30.3         lpSolve_5.6.19             
##   [7] scatterpie_0.2.1            ScaledMatrix_1.4.1         
##   [9] munsell_0.5.0               codetools_0.2-19           
##  [11] preprocessCore_1.58.0       withr_2.5.2                
##  [13] colorspace_2.0-3            GOSemSim_2.22.0            
##  [15] Biobase_2.56.0              limSolve_1.5.7             
##  [17] highr_0.10                  knitr_1.45                 
##  [19] rstudioapi_0.15.0           SingleCellExperiment_1.18.1
##  [21] stats4_4.2.0                ggsignif_0.6.4             
##  [23] DOSE_3.22.1                 MatrixGenerics_1.8.1       
##  [25] GenomeInfoDbData_1.2.8      KMsurv_0.1-5               
##  [27] polyclip_1.10-6             bit64_4.0.5                
##  [29] farver_2.1.1                rhdf5_2.40.0               
##  [31] downloader_0.4              vctrs_0.6.4                
##  [33] treeio_1.21.2               generics_0.1.3             
##  [35] xfun_0.40                   timechange_0.2.0           
##  [37] R6_2.5.1                    doParallel_1.0.17          
##  [39] GenomeInfoDb_1.34.4         clue_0.3-61                
##  [41] graphlayouts_1.0.1          rsvd_1.0.5                 
##  [43] locfit_1.5-9.8              rhdf5filters_1.8.0         
##  [45] bitops_1.0-7                cachem_1.0.6               
##  [47] fgsea_1.22.0                gridGraphics_0.5-1         
##  [49] DelayedArray_0.22.0         assertthat_0.2.1           
##  [51] scales_1.2.1                ggraph_2.1.0               
##  [53] enrichplot_1.16.2           gtable_0.3.4               
##  [55] beachmat_2.12.0             tidygraph_1.2.3            
##  [57] rlang_1.1.1                 genefilter_1.78.0          
##  [59] GlobalOptions_0.1.2         splines_4.2.0              
##  [61] rstatix_0.7.2               lazyeval_0.2.2             
##  [63] broom_1.0.5                 yaml_2.3.7                 
##  [65] reshape2_1.4.4              abind_1.4-5                
##  [67] backports_1.4.1             qvalue_2.28.0              
##  [69] tools_4.2.0                 bookdown_0.36              
##  [71] ggplotify_0.1.2             jquerylib_0.1.4            
##  [73] RColorBrewer_1.1-3          proxy_0.4-27               
##  [75] BiocGenerics_0.42.0         Rcpp_1.0.9                 
##  [77] plyr_1.8.7                  sparseMatrixStats_1.8.0    
##  [79] zlibbioc_1.42.0             purrr_1.0.2                
##  [81] RCurl_1.98-1.7              GetoptLong_1.0.5           
##  [83] viridis_0.6.4               cowplot_1.1.1              
##  [85] S4Vectors_0.34.0            zoo_1.8-10                 
##  [87] SummarizedExperiment_1.26.1 ggrepel_0.9.1              
##  [89] cluster_2.1.3               fs_1.5.2                   
##  [91] magrittr_2.0.3              data.table_1.14.2          
##  [93] DO.db_2.9                   circlize_0.4.15            
##  [95] matrixStats_0.62.0          GSVA_1.44.5                
##  [97] evaluate_0.22               xtable_1.8-4               
##  [99] XML_3.99-0.14               IRanges_2.30.0             
## [101] gridExtra_2.3               shape_1.4.6                
## [103] compiler_4.2.0              shadowtext_0.1.2           
## [105] crayon_1.5.2                htmltools_0.5.6.1          
## [107] ggfun_0.1.3                 tidyr_1.3.0                
## [109] geneplotter_1.74.0          aplot_0.2.2                
## [111] lubridate_1.9.3             DBI_1.1.3                  
## [113] tweenr_2.0.2                corrplot_0.92              
## [115] MASS_7.3-60                 Matrix_1.6-3               
## [117] car_3.1-2                   cli_3.6.1                  
## [119] quadprog_1.5-8              parallel_4.2.0             
## [121] igraph_1.3.2                GenomicRanges_1.48.0       
## [123] pkgconfig_2.0.3             km.ci_0.5-6                
## [125] foreach_1.5.2               ggtree_3.4.4               
## [127] annotate_1.74.0             bslib_0.5.1                
## [129] XVector_0.36.0              yulab.utils_0.1.0          
## [131] stringr_1.5.0               digest_0.6.29              
## [133] graph_1.74.0                Biostrings_2.64.0          
## [135] rmarkdown_2.25              fastmatch_1.1-4            
## [137] survMisc_0.5.6              tidytree_0.4.5             
## [139] dendextend_1.17.1           DelayedMatrixStats_1.18.2  
## [141] GSEABase_1.58.0             rjson_0.2.21               
## [143] lifecycle_1.0.3             nlme_3.1-157               
## [145] jsonlite_1.8.0              Rhdf5lib_1.18.2            
## [147] carData_3.0-5               viridisLite_0.4.2          
## [149] limma_3.52.4                fansi_1.0.3                
## [151] pillar_1.9.0                lattice_0.20-45            
## [153] KEGGREST_1.36.3             fastmap_1.1.1              
## [155] httr_1.4.7                  GO.db_3.15.0               
## [157] emo_0.0.0.9000              glue_1.6.2                 
## [159] png_0.1-7                   iterators_1.0.14           
## [161] glmnet_4.1-8                bit_4.0.5                  
## [163] HDF5Array_1.24.2            ggforce_0.4.1              
## [165] class_7.3-22                stringi_1.7.6              
## [167] sass_0.4.7                  blob_1.2.4                 
## [169] BiocSingular_1.12.0         DESeq2_1.36.0              
## [171] memoise_2.0.1               irlba_2.3.5                
## [173] tidyverse_2.0.0             e1071_1.7-13               
## [175] ape_5.6-2
```
