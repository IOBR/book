
# **How to use IOBR**

## 🔍 The main pipeline of IOBR

<div class="figure" style="text-align: center">
<img src="./fig/IOBR-Package.png" alt="The main pipeline of IOBR" width="95%" />
<p class="caption">(\#fig:flowchart)The main pipeline of IOBR</p>
</div>

## 👈 Main Functions of IOBR

* <div style="color:green">**Data Preparation: data annotation and transformation**</div> 
   * `count2tpm()`: transform gene expression count data into Transcripts Per Million (TPM) values. This function supports gene IDs of type "Ensembl", "Entrez", or "Symbol", and retrieves gene length information using either an online connection to the bioMart database or a local dataset (specified by the source parameter).
   * `anno_eset()`: annotate an ExpressionSet object (eset) with gene symbols using the provided annotation data.  It retains only the rows with probes that have matching identifiers in the annotation data. The function handles duplicates according to the specified method. The output is an annotated and cleaned expression set.
   * `remove_duplicate_genes()`: remove duplicate gene symbols from gene expression data. The retention of gene symbols is based either on their mean values (if method is set as "mean") or standard deviation values (if method is set as "sd").
   * `mouse2human_eset()`: convert mouse gene symbols to human gene symbols of expression set.
   * `find_outlier_samples()`: analyze gene expression data and identify potential outlier samples based on connectivity analysis. By utilizing the "WGCNA" package, this function calculates the normalized adjacency and connectivity z-scores for each sample. It also offers multiple parameters to customize analysis and visualization.
   * `remove_batcheffect()`: remove batch effects from given expression datasets and visualize the corrected data using principal component analysis (PCA). It takes three expression datasets as input and performs batch effect correction using the "sva::ComBat" or "sva::ComBat_seq" methods. The function then generates PCA plots to compare the data before and after correction.
   </br>

* <div style="color:green">**TME Deconvolution Module: integrate multiple algorithms to decode immune contexture**</div> 
   * `deconvo_tme()`: decode the TME infiltration using various deconvolution methodologies, based on bulk RNAseq, microarray or single cell RNAseq data. It currently supports methods include "CIBERSORT", "MCPcounter", "EPIC", "xCell", "IPS", "estimate", "quanTIseq", "TIMER", "SVR" and "lsei".
   * `generateRef()`: generate a novel gene reference data for specific feature deconvolution, such as infiltrating cell, utilizing different methods to identify differentially expressed genes (DEGs) . The function supports both "limma" and "DESeq2" methods.The resulting gene reference data can be used for `deconvo_tme()` with the "SVR" and "lsei" algorithms.
   * `generateRef_seurat()`: take a Seurat object "sce" and additional parameters to perform various operations for generating reference gene expression data. It allows for specifying cell types, proportions, assays, preprocessing options, and statistical testing parameters.The resulting gene reference data can be used for `deconvo_tme()` with the "svr" and "lsei" algorithms.
</br>

* <div style="color:green">**Signature Module: calculatint signature scores, estimate phenotype related signatures and corresponding genes, and evaluate signatures generated from single-cell RNA sequencing data **</div>
  * `calculate_sig_score()`: estimate the interested signatures enrolled in IOBR R package, which involves TME-associated, tumor-metabolism, and tumor-intrinsic signatures.The supported methods for signature score calculation include "PCA", "ssGSEA", "z-score", and Integration.
  * `feature_manipulation()`: manipulate features including the cell fraction and signatures generated from multi-omics data for latter analysis and model construction. Remove missing values, outliers and variables without significant variance.
  * `format_signatures()`: generate the object for `calculate_sig_score()` function, by inputting a data frame with signatures as column names of corresponding gene sets, and return a list contain the signature information for calculating multiple signature scores.
  * `format_msigdb()`: transform the signature gene sets data  with gmt format, which is not included in the signature collection and might be downloaded in the MSgiDB website, into the object of `calculate_sig_score()`function.
  * `sig_gsea()`: conduct Gene Set Enrichment Analysis (GSEA) to identify significant gene sets based on differential gene expression data. This function performs GSEA using the fgsea package and provides visualizations and results in the form of tables and plots. It supports the utilization of user-defined gene sets or the use of predefined gene sets from MSigDB. 
  * `get_sig_sc()`: get top gene signatures from single-cell differential analysis for `calculate_sig_score()` function. The input is a matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
</br>
    
* <div style="color:green">**Batch Analysis and Visualization: batch survival analysis and batch correlation analysis and other batch statistical analyses **</div>
  * `batch_surv()`: perform batch survival analysis. It calculates hazard ratios and confidence intervals for the specified variables based on the given data containing time-related information.
  * `subgroup_survival()`: extract hazard ratio and confidence intervals from a coxph object of subgroup analysis.
  * `batch_cor()`: batch analysis of correlation between two continuous variables using Pearson correlation coefficient or Spearman's rank correlation coefficient.
  * `batch_wilcoxon()`: perform Wilcoxon rank-sum tests on a given data set to compare the distribution of a specified feature between two groups. It computes the p-values and ranks the significant features based on the p-values. It returns a data frame with the feature names, p-values, adjusted p-values, logarithm of p-values, and a star rating based on the p-value ranges.
  * `batch_pcc()`: provide a batch way to calculate the partial correlation coefficient between feature and others when controlling a third variable.
  * `iobr_cor_plot()`: visualization of batch correlation analysis of signatures from "sig_group". Visualize the correlation between signature or phenotype with  expression of gene sets in target signature is also supported.
  * `cell_bar_plot()`: batch visualization of TME cell fraction, supporting input of deconvolution results from "CIBERSORT", "EPIC" and "quanTIseq" methodologies to further compare the TME cell distributions within one sample or among different samples.
  * `iobr_pca()`: perform Principal Component Analysis (PCA), which reduces the dimensionality of data while maintaining most of the original variance, and visualizes the PCA results on a scatter plot. 
  * `iobr_deg()`: perform differential expression analysis on gene expression data using the DESeq2 or limma method. It filters low count data, calculates fold changes and adjusted p-values, and identifies DEGs based on specified cutoffs. It also provides optional visualization tools such as volcano plots and heatmaps.
  * `get_cor()`: calculate and visualize the correlation between two variables in a dataset. It provides options to scale the data, handle missing values, and incorporate additional data. The function supports various correlation methods. It generates a correlation plot with optional subtypes or categories, including a regression line. 
  * `get_cor_matrix()`: calculate and visualize the correlation matrix between two sets of variables in a dataset. It provides flexibility in defining correlation methods, handling missing values, and incorporating additional data. The function supports various correlation methods, such as Pearson correlation, and displays the correlation result in a customizable plot.
  * `roc_time()`: generate a Receiver Operating Characteristic (ROC) plot over time to assess the predictive performance of one or more variables in survival analysis. It calculates the Area Under the Curve (AUC) for each specified time point and variable combination, and creates a multi-line ROC plot with corresponding AUC values annotated.
  * `sig_box()`: generate a boxplot with optional statistical comparisons. It takes in various parameters such as data, signature, variable, and more to customize the plot. It can be used to visualize and analyze data in a Seurat object or any other data frame.
  * `sig_heatmap()`: generate a heatmap plot based on input data, grouping variables, and optional conditions. The function allows customization of various parameters such as palette selection, scaling, color boxes, plot dimensions, and more. It provides flexibility in visualizing relationships between variables and groups in a concise and informative manner.
  * `sig_forest()`: create a forest plot for visualizing survival analysis results generated by "batch_surv".
  * `sig_roc()`: plot multiple ROC curves in a single graph, facilitating the comparison of different variables in terms of their ability to predict a binary response.
  * `sig_surv_plot()`: generat multiple Kaplan-Meier (KM) survival plots for a given signature or gene. It allows for detailed customization and is structured to handle various aspects of survival analysis.
  * `find_markers_in_bulk()`: find relevant results from the given gene expression data and meta information. It leverages the "Seurat" package to identify significant markers across multiple groups within the given data. The supported methods for comparison include "bootstrap", "delong" and  "venkatraman".
</br>

* <div style="color:green">**Signature Associated Mutation Module: identify and analyze mutations relevant to targeted signatures**</div>
  * `make_mut_matrix()`: transform the mutation data with MAF format(contain the columns of gene ID and the corresponding gene alterations which including SNP, indel and frameshift) into a mutation matrix in a suitable manner for further investigating signature relevant mutations.
  * `find_mutations()`: identify mutations associated with a distinct phenotype or signature. The function conducts the Cuzick test, Wilcoxon test, or both (when the method is set to "multi"). It generates box plots for the top genes identified through these statistical tests and creates oncoprints to graphically represent the mutation landscape across samples. 
</br>

* <div style="color:green">**Model Construction Module: feature selection and fast model construct to predict clinical phenotype**</div>
  * `BinomialModel()`: select features and construct a model to predict a binary phenotype. It accepts a dataset (x and y) as input and performs data processing, splitting into training and testing sets, and model fitting using both Lasso and Ridge regression techniques.
  * `PrognosticMode()`: select features and construct a model to predict clinical survival outcome. It primarily focuses on developing Lasso and Ridge regression models within the Cox proportional hazards framework.
  * `combine_pd_eset()`: combine the expression set (eset) with phenotype data (pdata).
  * `percent_bar_plot()`: create a percent bar plot based on the given data. The input is a data frame, with x and y-axis variables specified.
</br>


## 🗺 Current working environment


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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] rstudioapi_0.15.0 magrittr_2.0.3    knitr_1.45        timechange_0.2.0 
##  [5] R6_2.5.1          rlang_1.1.1       fastmap_1.1.1     highr_0.10       
##  [9] stringr_1.5.0     tools_4.2.0       xfun_0.40         cli_3.6.1        
## [13] jquerylib_0.1.4   htmltools_0.5.6.1 yaml_2.3.7        digest_0.6.29    
## [17] assertthat_0.2.1  lifecycle_1.0.3   crayon_1.5.2      bookdown_0.36    
## [21] purrr_1.0.2       sass_0.4.7        vctrs_0.6.4       glue_1.6.2       
## [25] cachem_1.0.6      evaluate_0.22     rmarkdown_2.25    emo_0.0.0.9000   
## [29] stringi_1.7.6     compiler_4.2.0    bslib_0.5.1       generics_0.1.3   
## [33] jsonlite_1.8.0    lubridate_1.9.3
```
