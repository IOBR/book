--- 
title: "IOBR (Immuno-Oncology Biological Research)"
author: "Dongqiang Zeng, Yiran Fang"
date: "2025-05-16"
output:
  pdf_document: default
  html_document:
    df_print: paged
documentclass: book
bibliography:
- book.bib
- packages.bib
url: https://iobr.github.io/book/
cover-image: ./fig/cover.png
geometry:
- top=1in
- left=1in
- right=1in
- bottom=1in
fontsize: 12pt
linestretch: 1.2
description: The tutorial of IOBR package
biblio-style: apalike
csl: "chicago-fullnote-bibliography.csl"
site: bookdown::bookdown_site
---
# **Introduction**{-}
Preface
<div class="figure" style="text-align: center">
<img src="./fig/IOBR-Workflow.png" alt="The workflow of IOBR" width="95%" />
<p class="caption">(\#fig:unnamed-chunk-1)The workflow of IOBR</p>
</div>

## 📚 Introduction

IOBR is the acronym for [Immuno-Oncology Biological Research](https://github.com/IOBR/IOBR).
Recent advances in next-generation sequencing have triggered the rapid accumulation of publicly available multi-omics data. The application of integrated omics to explore robust signatures for clinical translation is increasingly highlighted in immuno-oncology, but poses computational and biological challenges. This vignette aims to demonstrate how to use the package named IOBR to perform multi-omics immuno-oncology biological research to decipher tumour microenvironment and signatures for clinical translation.

This R package integrates 8 published methods for decoding the tumour microenvironment (TME) context: `CIBERSORT`, `TIMER`, `xCell`, `MCPcounter`, `ESITMATE`, `EPIC`, `IPS`, `quanTIseq`. In addition, 264 published signature gene sets have been collected by IOBR covering tumour microenvironment, tumour metabolism, m6A, exosomes, microsatellite instability and tertiary lymphoid structure. The `signature_collection_citation` function is run to obtain the source papers, and the `signature_collection` function returns the detailed signature genes of all given signatures. IOBR then uses three computational methods to calculate the signature score, including `PCA`, `z-score` and `ssGSEA`. Note that IOBR collected and used several approaches for variable transition, visualisation, batch survival analysis, feature selection and statistical analysis. Batch analysis and visualisation of results are supported. The details of how IOBR works are described below.


## 💿 License

**IOBR** was released under the GPL v3.0 license. See [LICENSE](https://github.com/IOBR/IOBR/blob/master/LICENSE) for details. The code contained in this book is simultaneously available under the [GPL license](https://www.gnu.org/licenses/why-not-lgpl.html); this means that you are free to use it in your packages, as long as you cite the source. The online version of this book is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.](https://creativecommons.org/licenses/by-nc-sa/4.0/)

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>


## 🏆 Citations

Zeng DQ, Fang YR, … , Yu GC, Liao WJ. Enhancing Immuno-Oncology Investigations Through Multidimensional Decoding of Tumour Microenvironment with IOBR 2.0, **Cell Reports Methods**, 2024 [https://doi.org/10.1016/j.crmeth.2024.100910](https://doi.org/10.1016/j.crmeth.2024.100910)

Fang YR, …, Liao WJ, Zeng DQ, Systematic Investigation of Tumor Microenvironment and Antitumor Immunity With IOBR, **Med Research**, 2025, [DOI:10.1002/mdr2.70001](https://onlinelibrary.wiley.com/doi/epdf/10.1002/mdr2.70001)

## ⏳ Major Updates

<div class="figure" style="text-align: center">
<img src="./fig/IOBR-Workflow2.png" alt="The workflow of IOBR" width="95%" />
<p class="caption">(\#fig:unnamed-chunk-2)The workflow of IOBR</p>
</div>

## ✉️ Reporting bugs
Please report bugs to the [Github issues page](https://github.com/IOBR/IOBR/issues)

E-mail any questions to  Dr. Fang <fyr_nate@163.com> or Dr. Zeng <interlaken@smu.edu.cn>



