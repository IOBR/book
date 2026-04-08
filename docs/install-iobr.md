
# **How to install IOBR**

## ⏳ Install Dependency Packages

It is essential that you have R 3.6.3 or above already installed on your computer or server (R>4.1 is recommended for some suggested packages). IOBR is a pipeline that utilizes many other R packages that are currently available from CRAN, Bioconductor (and GitHub).


``` r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
```

## 💻 Install IOBR package

When the dependent environments are built, users are able to install IOBR from github by typing the following code into your R session:

``` r
if (!requireNamespace("IOBR", quietly = TRUE)) pak::pak("IOBR/IOBR")

library(IOBR)
```


## 📌 How to update IOBR

``` r
pak::pak("IOBR/IOBR")
```


More see [IOBR README](https://github.com/IOBR/IOBR).
