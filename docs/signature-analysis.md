

# **Signature and relevant phenotypes**

## Loading packages
Load the IOBR package in your R session after the installation is complete:

```r
library(IOBR)
library(survminer)
library(tidyverse)
```

## Downloading data for example
Obtaining data set from GEO [Gastric cancer: GSE62254](https://pubmed.ncbi.nlm.nih.gov/25894828/) using `GEOquery` R package.

```r
if (!requireNamespace("GEOquery", quietly = TRUE))  BiocManager::install("GEOquery")
library("GEOquery")
# NOTE: This process may take a few minutes which depends on the internet connection speed. Please wait for its completion.
eset_geo <- getGEO(GEO = "GSE62254", getGPL  = F, destdir = "./")
eset    <-eset_geo[[1]]
eset    <-exprs(eset)
eset[1:5,1:5]
```

```
##           GSM1523727 GSM1523728 GSM1523729 GSM1523744 GSM1523745
## 1007_s_at  3.2176645  3.0624323  3.0279131   2.921683  2.8456013
## 1053_at    2.4050109  2.4394879  2.2442708   2.345916  2.4328582
## 117_at     1.4933412  1.8067380  1.5959665   1.839822  1.8326058
## 121_at     2.1965561  2.2812181  2.1865556   2.258599  2.1874363
## 1255_g_at  0.8698382  0.9502466  0.8125414   1.012860  0.9441993
```

Annotation of genes in the expression matrix and removal of duplicate genes.

```r
# Load the annotation file `anno_hug133plus2` in IOBR.
head(anno_hug133plus2)
```

```
## # A tibble: 6 × 2
##   probe_id  symbol 
##   <fct>     <fct>  
## 1 1007_s_at MIR4640
## 2 1053_at   RFC2   
## 3 117_at    HSPA6  
## 4 121_at    PAX8   
## 5 1255_g_at GUCA1A 
## 6 1294_at   MIR5193
```

```r
# Conduct gene annotation using `anno_hug133plus2` file; If identical gene symbols exists, these genes would be ordered by the mean expression levels. The gene symbol with highest mean expression level is selected and remove others. 

eset<-anno_eset(eset       = eset,
                annotation = anno_hug133plus2,
                symbol     = "symbol",
                probe      = "probe_id",
                method     = "mean")
eset[1:5, 1:3]
```

```
##              GSM1523727 GSM1523728 GSM1523729
## SH3KBP1        4.327974   4.316195   4.351425
## RPL41          4.246149   4.246808   4.257940
## EEF1A1         4.293762   4.291038   4.262199
## COX2           4.250288   4.283714   4.270508
## LOC101928826   4.219303   4.219670   4.213252
```


## Signature score estimation

### Signature collection of IOBR


```r
# Return available parameter options of signature estimation.
signature_score_calculation_methods
```

```
##           PCA        ssGSEA       z-score   Integration 
##         "pca"      "ssgsea"      "zscore" "integration"
```

```r
#TME associated signatures
names(signature_tme)[1:20]
```

```
##  [1] "CD_8_T_effector"            "DDR"                       
##  [3] "APM"                        "Immune_Checkpoint"         
##  [5] "CellCycle_Reg"              "Pan_F_TBRs"                
##  [7] "Histones"                   "EMT1"                      
##  [9] "EMT2"                       "EMT3"                      
## [11] "WNT_target"                 "FGFR3_related"             
## [13] "Cell_cycle"                 "Mismatch_Repair"           
## [15] "Homologous_recombination"   "Nucleotide_excision_repair"
## [17] "DNA_replication"            "Base_excision_repair"      
## [19] "TMEscoreA_CIR"              "TMEscoreB_CIR"
```


```r
#Metabolism related signatures
names(signature_metabolism)[1:20]
```

```
##  [1] "Cardiolipin_Metabolism"                    
##  [2] "Cardiolipin_Biosynthesis"                  
##  [3] "Cholesterol_Biosynthesis"                  
##  [4] "Citric_Acid_Cycle"                         
##  [5] "Cyclooxygenase_Arachidonic_Acid_Metabolism"
##  [6] "Prostaglandin_Biosynthesis"                
##  [7] "Purine_Biosynthesis"                       
##  [8] "Pyrimidine_Biosynthesis"                   
##  [9] "Dopamine_Biosynthesis"                     
## [10] "Epinephrine_Biosynthesis"                  
## [11] "Norepinephrine_Biosynthesis"               
## [12] "Fatty_Acid_Degradation"                    
## [13] "Fatty_Acid_Elongation"                     
## [14] "Fatty_Acid_Biosynthesis"                   
## [15] "Folate_One_Carbon_Metabolism"              
## [16] "Folate_biosynthesis"                       
## [17] "Gluconeogenesis"                           
## [18] "Glycolysis"                                
## [19] "Glycogen_Biosynthesis"                     
## [20] "Glycogen_Degradation"
```


```r
#Signatures associated with biomedical basic research: such as m6A and exosomes
names(signature_tumor)
```

```
##  [1] "Nature_metabolism_Hypoxia"                
##  [2] "Winter_hypoxia_signature"                 
##  [3] "Hu_hypoxia_signature"                     
##  [4] "Molecular_Cancer_m6A"                     
##  [5] "MT_exosome"                               
##  [6] "SR_exosome"                               
##  [7] "Positive_regulation_of_exosomal_secretion"
##  [8] "Negative_regulation_of_exosomal_secretion"
##  [9] "Exosomal_secretion"                       
## [10] "Exosome_assembly"                         
## [11] "Extracellular_vesicle_biogenesis"         
## [12] "MC_Review_Exosome1"                       
## [13] "MC_Review_Exosome2"                       
## [14] "CMLS_Review_Exosome"                      
## [15] "Ferroptosis"                              
## [16] "EV_Cell_2020"
```


```r
#signature collection including all aforementioned signatures 
names(signature_collection)[1:20]
```

```
##  [1] "CD_8_T_effector"            "DDR"                       
##  [3] "APM"                        "Immune_Checkpoint"         
##  [5] "CellCycle_Reg"              "Pan_F_TBRs"                
##  [7] "Histones"                   "EMT1"                      
##  [9] "EMT2"                       "EMT3"                      
## [11] "WNT_target"                 "FGFR3_related"             
## [13] "Cell_cycle"                 "Mismatch_Repair"           
## [15] "Homologous_recombination"   "Nucleotide_excision_repair"
## [17] "DNA_replication"            "Base_excision_repair"      
## [19] "TMEscoreA_CIR"              "TMEscoreB_CIR"
```



```r
#citation of signatures
signature_collection_citation[1:20, ]
```

```
## # A tibble: 20 × 6
##    Signatures                 `Published year` Journal         Title PMID  DOI  
##    <chr>                                 <dbl> <chr>           <chr> <chr> <chr>
##  1 CD_8_T_effector                        2018 Nature          TGFβ… 2944… 10.1…
##  2 DDR                                    2018 Nature          TGFβ… 2944… 10.1…
##  3 APM                                    2018 Nature          TGFβ… 2944… 10.1…
##  4 Immune_Checkpoint                      2018 Nature          TGFβ… 2944… 10.1…
##  5 CellCycle_Reg                          2018 Nature          TGFβ… 2944… 10.1…
##  6 Pan_F_TBRs                             2018 Nature          TGFβ… 2944… 10.1…
##  7 Histones                               2018 Nature          TGFβ… 2944… 10.1…
##  8 EMT1                                   2018 Nature          TGFβ… 2944… 10.1…
##  9 EMT2                                   2018 Nature          TGFβ… 2944… 10.1…
## 10 EMT3                                   2018 Nature          TGFβ… 2944… 10.1…
## 11 WNT_target                             2018 Nature          TGFβ… 2944… 10.1…
## 12 FGFR3_related                          2018 Nature          TGFβ… 2944… 10.1…
## 13 Cell_cycle                             2018 Nature          TGFβ… 2944… 10.1…
## 14 Mismatch_Repair                        2018 Nature          TGFβ… 2944… 10.1…
## 15 Homologous_recombination               2018 Nature          TGFβ… 2944… 10.1…
## 16 Nucleotide_excision_repair             2018 Nature          TGFβ… 2944… 10.1…
## 17 DNA_replication                        2018 Nature          TGFβ… 2944… 10.1…
## 18 Base_excision_repair                   2018 Nature          TGFβ… 2944… 10.1…
## 19 TMEscoreA_CIR                          2019 Cancer Immunol… Tumo… 3084… 10.1…
## 20 TMEscoreB_CIR                          2019 Cancer Immunol… Tumo… 3084… 10.1…
```

Three methodologies were adopted in the process of signature score evaluation, comprising Single-sample Gene Set Enrichment Analysis (ssGSEA), Principal component analysis (PCA), and Z-score.

### Estimated by PCA method

```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)

sig_tme <- t(column_to_rownames(sig_tme, var = "ID"))
sig_tme[1:5, 1:3]
```

```
##                   GSM1523727 GSM1523728 GSM1523729
## CD_8_T_effector   -2.5513794  0.7789141 -2.1770675
## DDR               -0.8747614  0.7425162 -1.3272054
## APM                1.1098368  2.1988688 -0.9516419
## Immune_Checkpoint -2.3701787  0.9455120 -1.4844104
## CellCycle_Reg      0.1063358  0.7583302 -0.3649795
```

### Estimated by ssGSEA methodology

This method is suitable for gene sets with a large number of genes, such as those of [GO, KEGG, REACTOME gene sets](https://www.gsea-msigdb.org/gsea/msigdb).

<div class="figure" style="text-align: center">
<img src="./fig/gsea.png" alt="Gene sets of MSigDb" width="95%" />
<p class="caption">(\#fig:unnamed-chunk-10)Gene sets of MSigDb</p>
</div>


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = go_bp,
                             method          = "ssgsea",
                             mini_gene_count = 2)
```

### Estimated by zscore function


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "zscore",
                             mini_gene_count = 2)
```

### Reference

**ssgsea**: Barbie, D.A. et al (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112.

**gsva**: Hänzelmann, S., Castelo, R. and Guinney, J. (2013). GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14(1):7.

**zscore**: Lee, E. et al (2008). Inferring pathway activity toward precise disease classification. PLoS Comp Biol, 4(11):e1000217.

**PCA method**: Mariathasan S, Turley SJ, Nickles D, et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548.


## Identifying features associated with survival


```r
data("pdata_acrg")
input <- combine_pd_eset(eset = sig_tme, pdata = pdata_acrg, scale = T)
res<- batch_surv(pdata    = input,
                 time     = "OS_time", 
                 status   = "OS_status", 
                 variable = colnames(input)[69:ncol(input)])
head(res)
```

```
## # A tibble: 6 × 5
##   ID                           P    HR CI_low_0.95 CI_up_0.95
##   <chr>                    <dbl> <dbl>       <dbl>      <dbl>
## 1 Folate_biosynthesis   1.00e-10 0.579       0.490      0.683
## 2 TMEscore_CIR          1.32e- 9 0.640       0.554      0.739
## 3 Glycogen_Biosynthesis 3.24e- 9 1.52        1.32       1.74 
## 4 Pan_F_TBRs            6.33e- 9 1.55        1.34       1.80 
## 5 TMEscoreB_CIR         7.17e- 9 1.52        1.32       1.75 
## 6 TMEscore_plus         8.08e- 9 0.638       0.547      0.743
```

```r
res<- res[nchar(res$ID)<=28, ]
p1<- sig_forest(res, signature = "ID", n = 20)
```

<img src="signature-analysis_files/figure-epub3/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

## Visulization using heatmap

Signatures和分子分型之间的关系
使用`IOBR`的`sig_heatmap`进行热图的可视化

```r
p2 <- sig_heatmap(input         = input, 
                  features      = res$ID[1:20],
                  group         = "Subtype", 
                  palette_group = "jama", 
                  palette       = 6)
```

<img src="signature-analysis_files/figure-epub3/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

## Focus on target signatures


```r
p1 <- sig_box(data           = input, 
              signature      = "Glycogen_Biosynthesis",
              variable       = "Subtype",
              jitter         = FALSE,
              cols           =  NULL,
              palette        = "jama",
              show_pvalue    = TRUE,
              size_of_pvalue = 5,
              hjust          = 1, 
              angle_x_text   = 60, 
              size_of_font   = 8)
```

```
## # A tibble: 6 × 8
##   .y.       group1    group2           p    p.adj p.format p.signif method  
##   <chr>     <chr>     <chr>        <dbl>    <dbl> <chr>    <chr>    <chr>   
## 1 signature EMT       MSI       5.39e-15 3.20e-14 5.4e-15  ****     Wilcoxon
## 2 signature EMT       MSS/TP53- 5.53e-13 2.8 e-12 5.5e-13  ****     Wilcoxon
## 3 signature EMT       MSS/TP53+ 1.90e-12 7.6 e-12 1.9e-12  ****     Wilcoxon
## 4 signature MSI       MSS/TP53- 1.14e- 3 3.4 e- 3 0.0011   **       Wilcoxon
## 5 signature MSI       MSS/TP53+ 7.05e- 3 1.4 e- 2 0.0071   **       Wilcoxon
## 6 signature MSS/TP53- MSS/TP53+ 7.16e- 1 7.2 e- 1 0.7161   ns       Wilcoxon
```

```r
p2 <- sig_box(data           = input, 
              signature      = "Pan_F_TBRs",
              variable       = "Subtype",
              jitter         = FALSE,
              cols           = NULL,
              palette        = "jama",
              show_pvalue    = TRUE,
              angle_x_text   = 60, 
              hjust          = 1, 
              size_of_pvalue = 5, 
              size_of_font   = 8)
```

```
## # A tibble: 6 × 8
##   .y.       group1    group2           p    p.adj p.format p.signif method  
##   <chr>     <chr>     <chr>        <dbl>    <dbl> <chr>    <chr>    <chr>   
## 1 signature EMT       MSI       7.98e-17 3.20e-16 <2e-16   ****     Wilcoxon
## 2 signature EMT       MSS/TP53- 1.70e-17 1   e-16 <2e-16   ****     Wilcoxon
## 3 signature EMT       MSS/TP53+ 2.57e-17 1.3 e-16 <2e-16   ****     Wilcoxon
## 4 signature MSI       MSS/TP53- 1.32e- 2 4   e- 2 0.013    *        Wilcoxon
## 5 signature MSI       MSS/TP53+ 6.99e- 2 1.4 e- 1 0.070    ns       Wilcoxon
## 6 signature MSS/TP53- MSS/TP53+ 4.02e- 1 4   e- 1 0.402    ns       Wilcoxon
```

```r
p3 <- sig_box(data           = input, 
              signature      = "Immune_Checkpoint",
              variable       = "Subtype",
              jitter          = TRUE,
              cols           = NULL,
              palette        = "jama",
              show_pvalue    = TRUE,
              angle_x_text   = 60, 
              hjust          = 1, 
              size_of_pvalue = 5, 
              size_of_font   = 8)
```

```
## # A tibble: 6 × 8
##   .y.       group1    group2           p        p.adj p.format p.signif method  
##   <chr>     <chr>     <chr>        <dbl>        <dbl> <chr>    <chr>    <chr>   
## 1 signature EMT       MSI       2.20e- 2 0.044        0.0220   *        Wilcoxon
## 2 signature EMT       MSS/TP53- 2.11e- 3 0.0085       0.0021   **       Wilcoxon
## 3 signature EMT       MSS/TP53+ 4.03e- 1 0.4          0.4026   ns       Wilcoxon
## 4 signature MSI       MSS/TP53- 9.13e-10 0.0000000055 9.1e-10  ****     Wilcoxon
## 5 signature MSI       MSS/TP53+ 5.03e- 4 0.0025       0.0005   ***      Wilcoxon
## 6 signature MSS/TP53- MSS/TP53+ 4.82e- 3 0.014        0.0048   **       Wilcoxon
```



```r
p1|p2|p3
```

<img src="signature-analysis_files/figure-epub3/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />


## Survival analysis
Signature的多种分层下的生存分析

```r
res <-       sig_surv_plot(input_pdata       = input, 
                           signature         = "Glycogen_Biosynthesis",
                           cols              = NULL, 
                           palette           = "jco",
                           project           = "ACRG",
                           time              = "OS_time",
                           status            = "OS_status",
                           time_type         = "month",
                           save_path         = "result")
```

```
##           ID   time status Glycogen_Biosynthesis group3 group2 bestcutoff
## 1 GSM1523727  88.73      0            -0.3612213 Middle    Low        Low
## 2 GSM1523728  88.23      0            -0.6926726    Low    Low        Low
## 3 GSM1523729  88.23      0            -0.9388531    Low    Low        Low
## 4 GSM1523744 105.70      0            -1.1825136    Low    Low        Low
## 5 GSM1523745 105.53      0            -0.3034304 Middle    Low        Low
## 6 GSM1523746  25.50      1             0.7517934   High   High       High
```

```
## [1] ">>>>>>>>>"
```

```r
res$plots
```

![](signature-analysis_files/figure-epub3/unnamed-chunk-17-1.png)<!-- -->


Signature在预测生存上的ROC

```r
p1<- roc_time(input      = input,  
             vars       = "Glycogen_Biosynthesis", 
             time       = "OS_time",
             status     = "OS_status", 
             time_point = c(12, 24, 36), 
             time_type  = "month",
             palette    = "jama",
             cols       = "normal",
             seed       = 1234, 
             show_col   = FALSE, 
             path       = "result", 
             main       = "OS",
             index      = 1,
             fig.type   = "pdf",
             width      = 5,
             height     = 5.2)
```

```
## [1] ">>>-- Range of Time: "
## [1]   1.0 105.7
```

```r
p2<- roc_time(input      = input,  
             vars       = "Glycogen_Biosynthesis", 
             time       = "RFS_time",
             status     = "RFS_status", 
             time_point = c(12, 24, 36), 
             time_type  = "month",
             palette    = "jama",
             cols       = "normal",
             seed       = 1234, 
             show_col   = FALSE, 
             path       = "result", 
             main       = "OS",
             index      = 1,
             fig.type   = "pdf",
             width      = 5,
             height     = 5.2)
```

```
## [1] ">>>-- Range of Time: "
## [1]   0.10 100.87
```

```r
p1|p2
```

![](signature-analysis_files/figure-epub3/unnamed-chunk-18-1.png)<!-- -->


## Batch correlation analysis 
寻找与目标signature相关的基因或者signatures

```r
res <- batch_cor(data = input, target = "Glycogen_Biosynthesis", feature = colnames(input)[69:ncol(input)])
head(res)
```

```
## # A tibble: 6 × 6
##   sig_names                         p.value statistic    p.adj log10pvalue stars
##   <chr>                               <dbl>     <dbl>    <dbl>       <dbl> <fct>
## 1 TMEscoreB_CIR                    8.89e-42     0.678 2.27e-39        41.1 **** 
## 2 Glycine__Serine_and_Threonine_M… 7.49e-40    -0.666 9.54e-38        39.1 **** 
## 3 Ether_Lipid_Metabolism           3.84e-39     0.662 3.27e-37        38.4 **** 
## 4 MDSC_Peng_et_al                  1.13e-38     0.659 7.21e-37        37.9 **** 
## 5 Glycerophospholipid_Metabolism   8.72e-38    -0.653 4.44e-36        37.1 **** 
## 6 TIP_Release_of_cancer_cell_anti… 2.32e-37    -0.650 9.86e-36        36.6 ****
```


```r
p1<- get_cor(eset = sig_tme, pdata = pdata_acrg, var1 = "Glycogen_Biosynthesis", var2 = "TMEscore_CIR", subtype = "Subtype", palette = "aaas")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  data[, var1] and data[, var2]
## S = 7282858, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.6184309 
## 
## [1] ">>>--- The exact p value is: 4.78971420439895e-33"
##       EMT       MSI MSS/TP53- MSS/TP53+ 
##        46        68       107        79
```

```r
p2<- get_cor(eset = sig_tme, pdata = pdata_acrg, var1 = "Glycogen_Biosynthesis", var2 = "TGFb.myCAF", subtype = "Subtype", palette = "aaas")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  data[, var1] and data[, var2]
## S = 2471758, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.4507143 
## 
## [1] ">>>--- The exact p value is: 2.04505761057615e-16"
##       EMT       MSI MSS/TP53- MSS/TP53+ 
##        46        68       107        79
```


```r
p1|p2
```

![](signature-analysis_files/figure-epub3/unnamed-chunk-21-1.png)<!-- -->



```r
feas1 <- c("Glycogen_Biosynthesis", "Ferroptosis")
feas2 <- c("Glutathione_Metabolism", "TMEscore_CIR", "Purine_Metabolism", "ICB_resistance_Peng_et_al", "Interleukins_Li_et_al", "TLS_Nature")
p <- get_cor_matrix(data           = input, 
                    feas1          = feas2, 
                    feas2          = feas1,
                    method         = "pearson",
                    font.size.star = 8, 
                    font.size      = 15, 
                    fill_by_cor    = FALSE, 
                    round.num      = 1)
```

![](signature-analysis_files/figure-epub3/unnamed-chunk-22-1.png)<!-- -->


## Visulization of correlations 

```r
input2 <- combine_pd_eset(eset = eset, pdata =  input[, c("ID", "Glycogen_Biosynthesis", "TLS_Nature", "Ferroptosis")])
feas1 <- c("Glycogen_Biosynthesis","TLS_Nature", "Ferroptosis")
feas2 <- signature_collection$CD_8_T_effector
feas2
```

```
## [1] "CD8A"   "GZMA"   "GZMB"   "IFNG"   "CXCL9"  "CXCL10" "PRF1"   "TBX21"
```

```r
p <- get_cor_matrix(data           = input2, 
                    feas1          = feas2, 
                    feas2          = feas1,
                    method         = "pearson",
                    scale          = T, 
                    font.size.star = 8, 
                    font.size      = 15, 
                    fill_by_cor    = FALSE, 
                    round.num      = 1)
```

![](signature-analysis_files/figure-epub3/unnamed-chunk-23-1.png)<!-- -->


```r
p <- get_cor_matrix(data           = input2, 
                    feas1          = feas2, 
                    feas2          = feas1,
                    method         = "pearson",
                    scale          = T, 
                    font.size.star = 8, 
                    font.size      = 15, 
                    fill_by_cor    = TRUE, 
                    round.num      = 2)
```

![](signature-analysis_files/figure-epub3/unnamed-chunk-24-1.png)<!-- -->

