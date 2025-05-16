
# **How to use IOBR**

## 🔍 The main pipeline of IOBR

<div class="figure" style="text-align: center">
<img src="./fig/IOBR-Package.png" alt="The main pipeline of IOBR" width="95%" />
<p class="caption">(\#fig:flowchart)The main pipeline of IOBR</p>
</div>

## 👆 Main Functions of IOBR

* <div style="color:green">**Data Preparation: data annotation and transformation**</div> 
   * **count2tpm()**: transform gene expression count data into Transcripts Per Million (TPM) values. This function supports gene IDs of type "Ensembl", "Entrez", or "Symbol", and retrieves gene length information using either an online connection to the bioMart database or a local dataset (specified by the source parameter).
   * **anno_eset()**: annotate an ExpressionSet object (eset) with gene symbols using the provided annotation data.  It retains only the rows with probes that have matching identifiers in the annotation data. The function handles duplicates according to the specified method. The output is an annotated and cleaned expression set.
   * **remove_duplicate_genes()**: remove duplicate gene symbols from gene expression data. The retention of gene symbols is based either on their mean values (if method is set as "mean") or standard deviation values (if method is set as "sd").
   * **mouse2human_eset()**: convert mouse gene symbols to human gene symbols of expression set.
   * **find_outlier_samples()**: analyze gene expression data and identify potential outlier samples based on connectivity analysis. By utilizing the "WGCNA" package, this function calculates the normalized adjacency and connectivity z-scores for each sample. It also offers multiple parameters to customize analysis and visualization.
   * **remove_batcheffect()**: remove batch effects from given expression datasets and visualize the corrected data using principal component analysis (PCA). It takes three expression datasets as input and performs batch effect correction using the "sva::ComBat" or "sva::ComBat_seq" methods. The function then generates PCA plots to compare the data before and after correction.
   </br>

* <div style="color:green">**TME Deconvolution Module: integrate multiple algorithms to decode immune contexture**</div> 
   * **deconvo_tme()**: decode the TME infiltration using various deconvolution methodologies, based on bulk RNAseq, microarray or single cell RNAseq data. It currently supports methods include "CIBERSORT", "MCPcounter", "EPIC", "xCell", "IPS", "estimate", "quanTIseq", "TIMER", "SVR" and "lsei".
   * **generateRef()**: generate a novel gene reference data for specific feature deconvolution, such as infiltrating cell, utilizing different methods to identify differentially expressed genes (DEGs) . The function supports both "limma" and "DESeq2" methods.The resulting gene reference data can be used for **deconvo_tme()** with the "svr" and "lsei" algorithms.
   * **generateRef_seurat()**: take a Seurat object "sce" and additional parameters to perform various operations for generating reference gene expression data. It allows for specifying cell types, proportions, assays, preprocessing options, and statistical testing parameters.The resulting gene reference data can be used for **deconvo_tme()** with the "svr" and "lsei" algorithms.
</br>

* <div style="color:green">**Signature Module: calculating signature scores, estimate phenotype related signatures and corresponding genes, and evaluate signatures generated from single-cell RNA sequencing data **</div>
  * **calculate_sig_score()**: estimate the interested signatures enrolled in IOBR R package, which involves TME-associated, tumor-metabolism, and tumor-intrinsic signatures.The supported methods for signature score calculation include "PCA", "ssGSEA", "z-score", and Integration.
  * **feature_manipulation()**: manipulate features including the cell fraction and signatures generated from multi-omics data for latter analysis and model construction. Remove missing values, outliers and variables without significant variance.
  * **format_signatures()**: generate the object for **calculate_sig_score()** function, by inputting a data frame with signatures as column names of corresponding gene sets, and return a list contain the signature information for calculating multiple signature scores.
  * **format_msigdb()**: transform the signature gene sets data  with gmt format, which is not included in the signature collection and might be downloaded in the MSgiDB website, into the object of **calculate_sig_score()** function.
  * **sig_gsea()**: conduct Gene Set Enrichment Analysis (GSEA) to identify significant gene sets based on differential gene expression data. This function performs GSEA using the fgsea package and provides visualizations and results in the form of tables and plots. It supports the utilization of user-defined gene sets or the use of predefined gene sets from MSigDB. 
  * **get_sig_sc()**: get top gene signatures from single-cell differential analysis for **calculate_sig_score()** function. The input is a matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
</br>
    
* <div style="color:green">**Batch Analysis and Visualization: batch survival analysis and batch correlation analysis and other batch statistical analyses **</div>
  * **batch_surv()**: perform batch survival analysis. It calculates hazard ratios and confidence intervals for the specified variables based on the given data containing time-related information.
  * **subgroup_survival()**: extract hazard ratio and confidence intervals from a coxph object of subgroup analysis.
  * **batch_cor()**: batch analysis of correlation between two continuous variables using Pearson correlation coefficient or Spearman's rank correlation coefficient.
  * **batch_wilcoxon()**: perform Wilcoxon rank-sum tests on a given data set to compare the distribution of a specified feature between two groups. It computes the p-values and ranks the significant features based on the p-values. It returns a data frame with the feature names, p-values, adjusted p-values, logarithm of p-values, and a star rating based on the p-value ranges.
  * **batch_pcc()**: provide a batch way to calculate the partial correlation coefficient between feature and others when controlling a third variable.
  * **iobr_cor_plot()**: visualization of batch correlation analysis of signatures from "sig_group". Visualize the correlation between signature or phenotype with  expression of gene sets in target signature is also supported.
  * **cell_bar_plot()**: batch visualization of TME cell fraction, supporting input of deconvolution results from "CIBERSORT", "EPIC" and "quanTIseq" methodologies to further compare the TME cell distributions within one sample or among different samples.
  * **iobr_pca()**: perform Principal Component Analysis (PCA), which reduces the dimensionality of data while maintaining most of the original variance, and visualizes the PCA results on a scatter plot. 
  * **iobr_deg()**: perform differential expression analysis on gene expression data using the DESeq2 or limma method. It filters low count data, calculates fold changes and adjusted p-values, and identifies DEGs based on specified cutoffs. It also provides optional visualization tools such as volcano plots and heatmaps.
  * **get_cor()**: calculate and visualize the correlation between two variables in a dataset. It provides options to scale the data, handle missing values, and incorporate additional data. The function supports various correlation methods. It generates a correlation plot with optional subtypes or categories, including a regression line. 
  * **get_cor_matrix()**: calculate and visualize the correlation matrix between two sets of variables in a dataset. It provides flexibility in defining correlation methods, handling missing values, and incorporating additional data. The function supports various correlation methods, such as Pearson correlation, and displays the correlation result in a customizable plot.
  * **roc_time()**: generate a Receiver Operating Characteristic (ROC) plot over time to assess the predictive performance of one or more variables in survival analysis. It calculates the Area Under the Curve (AUC) for each specified time point and variable combination, and creates a multi-line ROC plot with corresponding AUC values annotated.
  * **sig_box()**: generate a boxplot with optional statistical comparisons. It takes in various parameters such as data, signature, variable, and more to customize the plot. It can be used to visualize and analyze data in a Seurat object or any other data frame.
  * **sig_heatmap()**: generate a heatmap plot based on input data, grouping variables, and optional conditions. The function allows customization of various parameters such as palette selection, scaling, color boxes, plot dimensions, and more. It provides flexibility in visualizing relationships between variables and groups in a concise and informative manner.
  * **sig_forest()**: create a forest plot for visualizing survival analysis results generated by "batch_surv".
  * **sig_roc()**: plot multiple ROC curves in a single graph, facilitating the comparison of different variables in terms of their ability to predict a binary response.
  * **sig_surv_plot()**: generat multiple Kaplan-Meier (KM) survival plots for a given signature or gene. It allows for detailed customization and is structured to handle various aspects of survival analysis.
  * **find_markers_in_bulk()**: find relevant results from the given gene expression data and meta information. It leverages the "Seurat" package to identify significant markers across multiple groups within the given data. The supported methods for comparison include "bootstrap", "delong" and  "venkatraman".
</br>

* <div style="color:green">**Signature Associated Mutation Module: identify and analyze mutations relevant to targeted signatures**</div>
  * **make_mut_matrix()**: transform the mutation data with MAF format(contain the columns of gene ID and the corresponding gene alterations which including SNP, indel and frameshift) into a mutation matrix in a suitable manner for further investigating signature relevant mutations.
  * **find_mutations()**: identify mutations associated with a distinct phenotype or signature. The function conducts the Cuzick test, Wilcoxon test, or both (when the method is set to "multi"). It generates box plots for the top genes identified through these statistical tests and creates oncoprints to graphically represent the mutation landscape across samples. 
</br>

* <div style="color:green">**Model Construction Module: feature selection and fast model construct to predict clinical phenotype**</div>
  * **BinomialModel()**: select features and construct a model to predict a binary phenotype. It accepts a dataset (x and y) as input and performs data processing, splitting into training and testing sets, and model fitting using both Lasso and Ridge regression techniques.
  * **PrognosticMode()**: select features and construct a model to predict clinical survival outcome. It primarily focuses on developing Lasso and Ridge regression models within the Cox proportional hazards framework.
  * **combine_pd_eset()**: combine the expression set (eset) with phenotype data (pdata).
  * **percent_bar_plot()**: create a percent bar plot based on the given data. The input is a data frame, with x and y-axis variables specified.
</br>


## 🌐 Current working environment


``` r
library(IOBR)
```


``` r
sessionInfo()
```

```
## R version 4.4.1 (2024-06-14 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26100)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
## [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
## [3] LC_MONETARY=Chinese (Simplified)_China.utf8
## [4] LC_NUMERIC=C                               
## [5] LC_TIME=Chinese (Simplified)_China.utf8    
## 
## time zone: Asia/Shanghai
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] IOBR_0.99.0            GSVA_1.52.3            e1071_1.7-16          
##  [4] preprocessCore_1.66.0  limSolve_1.5.7.1       limma_3.60.6          
##  [7] lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1         
## [10] purrr_1.0.2            readr_2.1.5            tidyr_1.3.1           
## [13] tidyverse_2.0.0        survminer_0.4.9        patchwork_1.3.0       
## [16] clusterProfiler_4.12.6 tidyHeatmap_1.8.1      ComplexHeatmap_2.20.0 
## [19] survival_3.6-4         ggpubr_0.6.0           ggplot2_3.5.1         
## [22] dplyr_1.1.4            tibble_3.2.1          
## 
## loaded via a namespace (and not attached):
##   [1] SpatialExperiment_1.14.0    IRanges_2.38.1             
##   [3] R.methodsS3_1.8.2           GSEABase_1.66.0            
##   [5] goftest_1.2-3               HDF5Array_1.32.1           
##   [7] Biostrings_2.72.1           vctrs_0.6.5                
##   [9] spatstat.random_3.3-2       digest_0.6.37              
##  [11] png_0.1-8                   shape_1.4.6.1              
##  [13] proxy_0.4-27                ggrepel_0.9.6              
##  [15] deldir_2.0-4                parallelly_1.38.0          
##  [17] magick_2.8.5                MASS_7.3-60.2              
##  [19] reshape2_1.4.4              httpuv_1.6.15              
##  [21] foreach_1.5.2               BiocGenerics_0.50.0        
##  [23] qvalue_2.36.0               withr_3.0.2                
##  [25] xfun_0.48                   ggfun_0.1.7                
##  [27] memoise_2.0.1               gson_0.1.0                 
##  [29] tidytree_0.4.6              zoo_1.8-12                 
##  [31] GlobalOptions_0.1.2         pbapply_1.7-2              
##  [33] R.oo_1.26.0                 Formula_1.2-5              
##  [35] KEGGREST_1.44.1             promises_1.3.0             
##  [37] httr_1.4.7                  rstatix_0.7.2              
##  [39] rhdf5filters_1.16.0         globals_0.16.3             
##  [41] fitdistrplus_1.2-1          rhdf5_2.48.0               
##  [43] rstudioapi_0.17.1           UCSC.utils_1.0.0           
##  [45] miniUI_0.1.1.1              generics_0.1.3             
##  [47] DOSE_3.30.5                 S4Vectors_0.42.1           
##  [49] zlibbioc_1.50.0             ScaledMatrix_1.12.0        
##  [51] ggraph_2.2.1                polyclip_1.10-7            
##  [53] GenomeInfoDbData_1.2.12     quadprog_1.5-8             
##  [55] SparseArray_1.4.8           xtable_1.8-4               
##  [57] doParallel_1.0.17           evaluate_1.0.1             
##  [59] S4Arrays_1.4.1              hms_1.1.3                  
##  [61] glmnet_4.1-8                GenomicRanges_1.56.2       
##  [63] bookdown_0.41.1             irlba_2.3.5.1              
##  [65] colorspace_2.1-1            ROCR_1.0-11                
##  [67] reticulate_1.39.0           spatstat.data_3.1-2        
##  [69] magrittr_2.0.3              lmtest_0.9-40              
##  [71] later_1.3.2                 viridis_0.6.5              
##  [73] ggtree_3.12.0               lattice_0.22-6             
##  [75] spatstat.geom_3.3-3         future.apply_1.11.3        
##  [77] XML_3.99-0.17               scattermore_1.2            
##  [79] shadowtext_0.1.4            cowplot_1.1.3              
##  [81] matrixStats_1.4.1           RcppAnnoy_0.0.22           
##  [83] class_7.3-22                pillar_1.9.0               
##  [85] nlme_3.1-166                iterators_1.0.14           
##  [87] beachmat_2.20.0             compiler_4.4.1             
##  [89] RSpectra_0.16-2             stringi_1.8.4              
##  [91] tensor_1.5                  SummarizedExperiment_1.34.0
##  [93] dendextend_1.18.1           plyr_1.8.9                 
##  [95] crayon_1.5.3                abind_1.4-8                
##  [97] gridGraphics_0.5-1          locfit_1.5-9.10            
##  [99] sp_2.1-4                    graphlayouts_1.2.0         
## [101] bit_4.5.0                   fastmatch_1.1-4            
## [103] codetools_0.2-20            BiocSingular_1.20.0        
## [105] pec_2023.04.12              bslib_0.8.0                
## [107] GetoptLong_1.0.5            plotly_4.10.4              
## [109] mime_0.12                   splines_4.4.1              
## [111] circlize_0.4.16             Rcpp_1.0.13                
## [113] fastDummies_1.7.4           sparseMatrixStats_1.16.0   
## [115] knitr_1.48                  blob_1.2.4                 
## [117] utf8_1.2.4                  clue_0.3-65                
## [119] fs_1.6.4                    listenv_0.9.1              
## [121] ggsignif_0.6.4              ggplotify_0.1.2            
## [123] Matrix_1.7-0                statmod_1.5.0              
## [125] tzdb_0.4.0                  lpSolve_5.6.21             
## [127] tweenr_2.0.3                pkgconfig_2.0.3            
## [129] tools_4.4.1                 cachem_1.1.0               
## [131] RSQLite_2.3.7               viridisLite_0.4.2          
## [133] DBI_1.2.3                   numDeriv_2016.8-1.1        
## [135] fastmap_1.2.0               rmarkdown_2.28             
## [137] scales_1.3.0                ica_1.0-3                  
## [139] Seurat_5.1.0                broom_1.0.7                
## [141] sass_0.4.9                  dotCall64_1.2              
## [143] graph_1.82.0                carData_3.0-5              
## [145] RANN_2.6.2                  farver_2.1.2               
## [147] tidygraph_1.3.1             scatterpie_0.2.4           
## [149] yaml_2.3.10                 MatrixGenerics_1.16.0      
## [151] cli_3.6.3                   stats4_4.4.1               
## [153] leiden_0.4.3.1              lifecycle_1.0.4            
## [155] uwot_0.2.2                  Biobase_2.64.0             
## [157] mvtnorm_1.3-1               lava_1.8.0                 
## [159] backports_1.5.0             annotate_1.82.0            
## [161] BiocParallel_1.38.0         timechange_0.3.0           
## [163] gtable_0.3.6                rjson_0.2.23               
## [165] ggridges_0.5.6              progressr_0.14.0           
## [167] parallel_4.4.1              ape_5.8                    
## [169] jsonlite_1.8.9              RcppHNSW_0.6.0             
## [171] bit64_4.5.2                 assertthat_0.2.1           
## [173] Rtsne_0.17                  yulab.utils_0.1.7          
## [175] spatstat.utils_3.1-0        SeuratObject_5.0.2         
## [177] jquerylib_0.1.4             highr_0.11                 
## [179] GOSemSim_2.30.2             survMisc_0.5.6             
## [181] spatstat.univar_3.0-1       R.utils_2.12.3             
## [183] lazyeval_0.2.2              shiny_1.9.1                
## [185] htmltools_0.5.8.1           enrichplot_1.24.4          
## [187] KMsurv_0.1-5                GO.db_3.19.1               
## [189] sctransform_0.4.1           rappdirs_0.3.3             
## [191] glue_1.8.0                  timereg_2.0.6              
## [193] spam_2.11-0                 httr2_1.0.5                
## [195] XVector_0.44.0              treeio_1.28.0              
## [197] gridExtra_2.3               igraph_2.1.1               
## [199] R6_2.5.1                    SingleCellExperiment_1.26.0
## [201] DESeq2_1.44.0               km.ci_0.5-6                
## [203] cluster_2.1.6               Rhdf5lib_1.26.0            
## [205] aplot_0.2.3                 GenomeInfoDb_1.40.1        
## [207] DelayedArray_0.30.1         tidyselect_1.2.1           
## [209] ggforce_0.4.2               car_3.1-3                  
## [211] AnnotationDbi_1.66.0        future_1.34.0              
## [213] emo_0.0.0.9000              rsvd_1.0.5                 
## [215] munsell_0.5.1               KernSmooth_2.23-24         
## [217] data.table_1.16.2           htmlwidgets_1.6.4          
## [219] fgsea_1.30.0                RColorBrewer_1.1-3         
## [221] rlang_1.1.4                 spatstat.sparse_3.1-0      
## [223] spatstat.explore_3.3-3      fansi_1.0.6                
## [225] timeROC_0.4                 prodlim_2024.06.25
```

