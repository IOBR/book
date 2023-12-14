

# **Signature Score Calculation**

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

Signatures associated with basic biomedical research, such as m6A, TLS, ferroptosis and exosomes.

```r
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

`signature_collection` including all aforementioned signatures 

```r
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

The evaluation of signature scores involved three methodologies: Single-sample Gene Set Enrichment Analysis (ssGSEA), Principal Component Analysis (PCA), and Z-score.

### Estimation of signature using PCA method
The PCA method is ideal for gene sets with co-expression. Heatmaps and correlation matrices can be used to determine if co-expression is present in the applicable gene set.

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

### Estimated using the ssGSEA methodology

This method is appropriate for gene sets that contain a large number of genes (> 30 genes), such as those of [GO, KEGG, REACTOME gene sets](https://www.gsea-msigdb.org/gsea/msigdb).

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

### Calculated using the z-score function.


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "zscore",
                             mini_gene_count = 2)
```

### Calculated using all three methods at the same time


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "integration",
                             mini_gene_count = 2)
```
The same SIGNATURE in this case will be scored using all three methods simultaneously.

```r
colnames(sig_tme)[grep(colnames(sig_tme), pattern = "CD_8_T_effector")]
```

The `select_method()` function allows the user to extract data using various methods.

```r
sig_tme_pca <- select_method(data = sig_tme, method = "pca")
colnames(sig_tme_pca)[grep(colnames(sig_tme_pca), pattern = "CD_8_T_effector")]
```


### How to customise the signature gene list



### References

**ssgsea**: Barbie, D.A. et al (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112.

**gsva**: Hänzelmann, S., Castelo, R. and Guinney, J. (2013). GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14(1):7.

**zscore**: Lee, E. et al (2008). Inferring pathway activity toward precise disease classification. PLoS Comp Biol, 4(11):e1000217.

**PCA method**: Mariathasan S, Turley SJ, Nickles D, et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548.

**MSigDB**:Dolgalev I (2022). msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format. R  package version 7.5.1. (https://www.gsea-msigdb.org/gsea/msigdb/)

