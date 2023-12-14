
# **How to install IOBR**

## ⏳ Install Dependency Packages

It is essential that you have R 3.6.3 or above already installed on your computer or server. IOBR is a pipeline that utilizes many other R packages that are currently available from CRAN, Bioconductor and GitHub.


```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
          "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
          "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr', "PMCMRplus")
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
