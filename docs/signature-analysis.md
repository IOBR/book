

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
## # A tibble: 6 x 2
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
## # A tibble: 20 x 6
##    Signatures                 `Published year` Journal         Title PMID  DOI  
##    <chr>                                 <dbl> <chr>           <chr> <chr> <chr>
##  1 CD_8_T_effector                        2018 Nature          TGFβ~ 2944~ 10.1~
##  2 DDR                                    2018 Nature          TGFβ~ 2944~ 10.1~
##  3 APM                                    2018 Nature          TGFβ~ 2944~ 10.1~
##  4 Immune_Checkpoint                      2018 Nature          TGFβ~ 2944~ 10.1~
##  5 CellCycle_Reg                          2018 Nature          TGFβ~ 2944~ 10.1~
##  6 Pan_F_TBRs                             2018 Nature          TGFβ~ 2944~ 10.1~
##  7 Histones                               2018 Nature          TGFβ~ 2944~ 10.1~
##  8 EMT1                                   2018 Nature          TGFβ~ 2944~ 10.1~
##  9 EMT2                                   2018 Nature          TGFβ~ 2944~ 10.1~
## 10 EMT3                                   2018 Nature          TGFβ~ 2944~ 10.1~
## 11 WNT_target                             2018 Nature          TGFβ~ 2944~ 10.1~
## 12 FGFR3_related                          2018 Nature          TGFβ~ 2944~ 10.1~
## 13 Cell_cycle                             2018 Nature          TGFβ~ 2944~ 10.1~
## 14 Mismatch_Repair                        2018 Nature          TGFβ~ 2944~ 10.1~
## 15 Homologous_recombination               2018 Nature          TGFβ~ 2944~ 10.1~
## 16 Nucleotide_excision_repair             2018 Nature          TGFβ~ 2944~ 10.1~
## 17 DNA_replication                        2018 Nature          TGFβ~ 2944~ 10.1~
## 18 Base_excision_repair                   2018 Nature          TGFβ~ 2944~ 10.1~
## 19 TMEscoreA_CIR                          2019 Cancer Immunol~ Tumo~ 3084~ 10.1~
## 20 TMEscoreB_CIR                          2019 Cancer Immunol~ Tumo~ 3084~ 10.1~
```

The evaluation of signature scores involved three methodologies: Single-sample Gene Set Enrichment Analysis (ssGSEA), Principal Component Analysis (PCA), and Z-score.

## Estimation of signature using PCA method
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

## Estimated using the ssGSEA methodology

This method is appropriate for gene sets that contain a large number of genes (> 30 genes), such as those of [GO, KEGG, REACTOME gene sets](https://www.gsea-msigdb.org/gsea/msigdb).

\begin{figure}

{\centering \includegraphics[width=0.95\linewidth]{./fig/gsea} 

}

\caption{Gene sets of MSigDb}(\#fig:unnamed-chunk-10)
\end{figure}


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = go_bp,
                             method          = "ssgsea",
                             mini_gene_count = 2)
```

## Calculated using the z-score function.


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "zscore",
                             mini_gene_count = 2)
```

## Calculated using all three methods at the same time


```r
sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = eset,
                             signature       = signature_collection,
                             method          = "integration",
                             mini_gene_count = 2)
```
The same signature in this case will be scored using all three methods simultaneously.

```r
colnames(sig_tme)[grep(colnames(sig_tme), pattern = "CD_8_T_effector")]
```

The `select_method()` function allows the user to extract data using various methods.

```r
sig_tme_pca <- select_method(data = sig_tme, method = "pca")
colnames(sig_tme_pca)[grep(colnames(sig_tme_pca), pattern = "CD_8_T_effector")]
```

## Classification of signatures
As more signatures related to the tumour microenvironment were collected in IOBR, and we may continue to add gene signatures related to the tumour microenvironment in the future, we have made a basic classification of these signatures by combining them with our analysis experience. Users can compare the signatures in the same group during the analysis process to improve the reliability and consistency of the conclusions.


```r
sig_group[1:8]
```

## How to customise the signature gene list for `calculate_signature_score`

### Method-1: Use excel for storage and construction
Users can collect gene signatures using either an `Excel` or `CSV` file. The format should have the name of the signature in the first row, followed by the genes contained in each signature from the second row onwards. Once imported, the function `format_signature` can be used to transform the data into a gene list of signatures required for `calculate_signature_score`. To import the file into R, users can use the functions read.csv or read_excel. It is important to note here that the user needs to use the longest signature as a criterion and then replace all the vacant grids using NA, otherwise an error may be reported when reading into R.

Here we provide a sample data `sig_excel`, please refer to this format to construct the required csv or excel files.

```r
data("sig_excel", package = "IOBR")
sig <- format_signatures(sig_excel)
print(sig[1:5])
```

```
## $Tcell_co_inhibitors
##  [1] "ADORA2A"  "BTLA"     "BTN2A2"   "BTN3A1"   "BTN3A2"   "BTNL2"   
##  [7] "C10orf54" "CSF1R"    "HAVCR2"   "IDO1"     "IL10"     "IL10RB"  
## [13] "KDR"      "KIR2DL1"  "SLAMF7"   "TGFB1"    "TIGIT"    "VRCN1"   
## [19] "VTCN1"    "CD247"    "CTLA4"    "CD160"    "CD244"    "CD274"   
## [25] "CD276"    "CD48"     "CD96"     "KIR2DL2"  "KIR2DL3"  "LAG3"    
## [31] "LAIR1"    "LGALS9"   "PVRL2"    "PDCD1"    "PDCD1LG2"
## 
## $Tcell_co_stimuiations
##  [1] "BTNL8"     "CD226"     "CD27"      "CD28"      "CD40"      "CD58"     
##  [7] "CD70"      "SLAMF1"    "TMIGD2"    "TNFRSF13B" "TNFRSF13C" "TNFRSF14" 
## [13] "TNFRSF4"   "TNFRSF8"   "TNFSF8"    "TNFSF9"    "ENTPD1"    "NT5E"     
## [19] "ICOS"      "TNFSF4"    "TNFSF15"   "CD80"      "CD86"      "EGFR"     
## [25] "HAVCR1"    "TNFSF18"   "ICOSLG"    "TNFSF13B"  "TNFRSF9"   "TNFSF13"  
## 
## $Tcell_function
## [1] "CD3E"  "CD4"   "CD8B"  "FOXP3" "GZMB"  "PRF1"  "TBX21" "IL2RA" "IKZF2"
## 
## $Tcell_checkpoint
##  [1] "CD274"    "CTLA4"    "LAG3"     "TIM3"     "TNFRSF9"  "TIGIT"   
##  [7] "CD226"    "CD7"      "GZMB"     "PRF1"     "TNFRSF18" "TNFRSF4" 
## [13] "HAVCR2"   "NLG1"     "CD4"      "CD8A"     "CD8B"     "FOXP3"   
## [19] "IL2"      "CXCL8"    "PDCD1"    "IFNG"    
## 
## $Teffctore_score
## [1] "CD8A"   "CXCL10" "CXCL9"  "GZMA"   "GZMB"   "IFNG"   "PRF1"   "TBX21"
```
For simple structures or when the number of signatures to be added is relatively small, the following two methods can also be used.

### Method-2: Build the list structure directly


```r
sig <- list("CD8" = c("CD8A",  "CXCL10", "CXCL9",  "GZMA",   "GZMB",   "IFNG",   "PRF1",   "TBX21"),
            "ICB" = c("CD274",   "PDCD1LG2", "CTLA4",    "PDCD1",    "LAG3",     "HAVCR2",   "TIGIT" ))
sig
```

```
## $CD8
## [1] "CD8A"   "CXCL10" "CXCL9"  "GZMA"   "GZMB"   "IFNG"   "PRF1"   "TBX21" 
## 
## $ICB
## [1] "CD274"    "PDCD1LG2" "CTLA4"    "PDCD1"    "LAG3"     "HAVCR2"   "TIGIT"
```
### Method3: Add the new signature to the existing gene list


```r
sig<- signature_tumor
sig$CD8 <- c("CD8A",  "CXCL10", "CXCL9",  "GZMA",   "GZMB",   "IFNG",   "PRF1",   "TBX21")
sig
```

```
## $Nature_metabolism_Hypoxia
##  [1] "ACOT7"  "SLC2A1" "ALDOA"  "CDKN3"  "ENO1"   "LDHA"   "MIF"    "MRPS17"
##  [9] "NDRG1"  "P4HA1"  "PGAM1"  "TPI1"   "TUBB6"  "VEGFA"  "ADM"   
## 
## $Winter_hypoxia_signature
## [1] "VEGF"  "GLUT1" "PDK-1" "EN01"  "HK2"   "CA9"   "AK3"   "CCNG2" "PFKB3"
## 
## $Hu_hypoxia_signature
##  [1] "FABP5"     "UCHL1"     "GAL"       "PLODDDIT4" "VEGF"      "ADM"      
##  [7] "ANGPTL4"   "NDRG1"     "NP"        "SLC16A3"   "C14ORF58"  "RRAGD"    
## 
## $Molecular_Cancer_m6A
##  [1] "METTL3"    "METTL14"   "RBM15"     "RBM15B"    "WTAP"      "KIAA1429" 
##  [7] "CBLL1"     "ZC3H13"    "ALKBH5"    "FTO"       "YTHDC1"    "YTHDC2"   
## [13] "YTHDF1"    "YTHDF2"    "YTHDF3"    "IGF2BP1"   "HNRNPA2B1" "HNRNPC"   
## [19] "FMR1"      "LRPPRC"    "ELAVL1"   
## 
## $MT_exosome
##  [1] "YWHAG"  "YWHAQ"  "CLTC"   "NCKAP1" "CFL1"   "ACTB"   "CCT4"   "RDX"   
##  [9] "GNA13"  "CTNNB1"
## 
## $SR_exosome
## [1] "HSP70" "HSP90" "CD9"   "CD63"  "CD81"  "CD82" 
## 
## $Positive_regulation_of_exosomal_secretion
##  [1] "ATP13A2" "CHMP2A"  "HGS"     "MYO5B"   "PDCD6IP" "RAB7"    "SDC1"   
##  [8] "SDC4"    "SDCBP"   "SMPD3"   "SNF8"    "STAM"    "TSG101"  "VPS4A"  
## 
## $Negative_regulation_of_exosomal_secretion
## [1] "VPS4B" "PRKN"  "RAB7" 
## 
## $Exosomal_secretion
## [1] "STEAP3" "TSG101" "RAB11A" "RAB27A" "COPS5" 
## 
## $Exosome_assembly
## [1] "CD34"    "PDCD6IP" "SDC1"    "SDC4"    "SDCBP"   "STAM"    "TSG101" 
## 
## $Extracellular_vesicle_biogenesis
##  [1] "ARRDC1"  "ARRDC4"  "ATP13A2" "CD34"    "CHMP2A"  "COPS5"   "HGS"    
##  [8] "MYO5B"   "PDCD6IP" "PRKN"    "RAB7"    "RAB11A"  "RAB27A"  "SDC1"   
## [15] "SDC4"    "SDCBP"   "SMPD3"   "SNF8"    "STAM"    "STEAP3"  "TSG101" 
## [22] "VPS4B"  
## 
## $MC_Review_Exosome1
##  [1] "TSG101"  "CD9"     "CD81"    "CD63"    "FLOT1"   "ITGB1"   "ITGA1"  
##  [8] "HSP70"   "AIP1"    "ALIX"    "PDCD6IP"
## 
## $MC_Review_Exosome2
##  [1] "RAB27A"  "RAB27B"  "PIKFYVE" "HRS"     "SYT7"    "CTTN"    "STAT3"  
##  [8] "PKM2"    "UNC13D"  "miR-155" "EGFR"    "RAS"     "EIF3C"   "LKB1"   
## [15] "STK11"  
## 
## $CMLS_Review_Exosome
##  [1] "HRS"      "STAM1"    "TSG101"   "CHMP4C"   "ALIX"     "VAT1"    
##  [7] "VPS4"     "CD9"      "CD82"     "CD63"     "LMP1"     "TSPAN8"  
## [13] "VAMP7"    "YKT6"     "PKM2"     "SNAP-23"  "RALA"     "RALB"    
## [19] "RAB2B"    "RAB5A"    "RAB9A"    "RAB7"     "RAB11"    "RAB27A"  
## [25] "RAB27B"   "RAB35"    "DGKA"     "PLD2"     "ARF6"     "ATG12"   
## [31] "ATG7"     "PIKFYVE"  "BST2"     "ATP6V0A4"
## 
## $Ferroptosis
##  [1] "ACSL4"      "AKR1C1-3"   "ALOXs"      "ATP5G3"     "CARS"      
##  [6] "CBS"        "CD44v"      "CHAC1"      "CISD1"      "CS"        
## [11] "DPP4"       "FANCD2"     "GCLC/GCLM"  "GLS2"       "GPX4"      
## [16] "GSS"        "HMGCR"      "HSPB1/5"    "KOD"        "LPCAT3"    
## [21] "MT1G"       "NCOA4"      "NFE2L2"     "PTGS2"      "RPL8"      
## [26] "SAT1"       "SLC7A11"    "SQS"        "TFRC"       "TP53"      
## [31] "TTC35/EMC2" "MESH1"     
## 
## $EV_Cell_2020
##  [1] "HSP90AB1" "HSP90AA1" "CD9"      "ALIX"     "FLOT1"    "FLOT2"   
##  [7] "TSG101"   "HSPA8"    "CD81"     "CD63"     "HBB"      "JCHAIN"  
## [13] "A2M"      "B2M"      "FN1"      "RAP1B"    "LGALS3BP" "GSN"     
## [19] "MSN"      "FLNA"     "ACTB"     "STOM"     "PRDX2"   
## 
## $CD8
## [1] "CD8A"   "CXCL10" "CXCL9"  "GZMA"   "GZMB"   "IFNG"   "PRF1"   "TBX21"
```

## How to export gene signature

Using the `output_sig` function, user can export the signatures of the list structure to a csv file for other purposes. This step is exactly the reverse of `format_signatures`.

```r
sig <- output_sig(signatures = signature_sc, format = "csv", file.name = "sc_signature")
sig[1:8, 1:5]
```

```
##   CD4_c0_Tcm CD4_c1_Treg CD4_c10_Tn_LEF1_ANKRD55 CD4_c11_Tisg CD4_c2_Tn
## 1      ANXA1       FOXP3                 ANKRD55        ISG15    NBEAL1
## 2       LMNA       IL2RA                    LEF1         IFI6      CCR7
## 3        VIM     TNFRSF4                    TCF7       IFI44L   GLTSCR2
## 4      KLRB1       TIGIT                   NOSIP          MX1      TCF7
## 5       IL7R      CARD16                    SELL        IFIT3    GNB2L1
## 6      ZFP36    TNFRSF18                   IL6ST        IFIT1      SELL
## 7    ZFP36L2        BATF                 LDLRAP1        RSAD2   C6orf48
## 8     GPR183       CTLA4                  RIPOR2        STAT1    TMEM66
```

## References

**ssgsea**: Barbie, D.A. et al (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112.

**gsva**: Hänzelmann, S., Castelo, R. and Guinney, J. (2013). GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14(1):7.

**zscore**: Lee, E. et al (2008). Inferring pathway activity toward precise disease classification. PLoS Comp Biol, 4(11):e1000217.

**PCA method**: Mariathasan S, Turley SJ, Nickles D, et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature. 2018 Feb 22;554(7693):544-548.

**MSigDB**:Dolgalev I (2022). msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format. R  package version 7.5.1. (https://www.gsea-msigdb.org/gsea/msigdb/)

