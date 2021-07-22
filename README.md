# Trying to plot Kaplan–Meier estimator for a given Ensembl Gene ID (ENSG) using Subjects with expression+clinical data in MMRF-COMMPASS

## install dependencies
``` r
install.packages("devtools")
devtools::install_github(c("RDocTaskForce/testextra","RDocTaskForce/parsetools","halpo/purrrogress"))
install.packages(c("BiocManager","shiny","markdown","knitr","sweave","xtable","DT","httpuv","sourcetools"))
BiocManager::install(c("EDASeq","genefilter","sva","limma","GenomicFeatures","GenomeInfoDb","GenomicRanges","SummarizedExperiment","EnsDb.Hsapiens.v79","S4Vectors","biomaRt","BiocStyle","edgeR","IRanges","TCGAbiolinks"))
``` 

## finally install [MMRFBiolinks](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab050/6209690)
``` 
devtools::install_github("marziasettino/MMRFBiolinks")
``` 

## download and process expression+clinical data in MMRF-COMMPASS
``` r
library(shiny)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(ggplot2)
library(TCGAbiolinks)
library(MMRFBiolinks)
MMRFclin <- MMRFGDC_QueryClinic(type = "clinical")
listSamples <- MMRFclin$bcr_patient_barcode
query <- GDCquery(project = "MMRF-COMMPASS",data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",experimental.strategy = "RNA-Seq",workflow.type="HTSeq - FPKM",barcode = listSamples)
GDCdownload(query, method = "api", files.per.chunk = 10)
MMRnaseqSE <- MMRFGDC_prepare(query,save = TRUE ,save.filename = "GDCdata.rda",directory = "GDCdata",summarizedExperiment = TRUE)
MMRFdataPrepro <- TCGAanalyze_Preprocessing(MMRnaseqSE)
save(file="MMRFdataPrepro.rds",MMRFdataPrepro)
save(file="MMRFclin.rds",MMRFclin)
savehistory("R.history")
``` 

## run shiny/app and provide [ENSG](https://www.ensembl.org/index.html) ID in GUI for gene/protein of interest, e.g. [ENSG00000196976](https://www.ensembl.org/Multi/Search/Results?q=ENSG00000196976;site=ensembl)
``` r
runApp('plotKM')
``` 

<!-- README.md is generated from README.Rmd. Please edit that file -->

# MMRFBiolinks

## An R package that extends TCGABiolink package for integrative analysis with MMRF-COMMPASS data

<!-- badges: start -->

<!-- badges: end -->

MMRFBiolinks extends TCGABiolink package for searching, downloading and
analyzing MMRF-COMMPASS data available at the NCI’s Genomic Data Commons
(GDC) Data Portal.

## Installation

Once R (version “4.0”) has been started, you can install the released
version of MMRFBiolinks from GitHub with:

``` r
devtools::install_github("marziasettino/MMRFBiolinks", build_vignettes = TRUE)
library(MMRFBiolinks)
```

## Required libraries

``` r
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(ggplot2)
```

## Vignettes

A list of all currently integrated vignettes can be obtained through:

``` r
vignette(package="MMRFBiolinks")
```

The best way to view vignettes is in your web browser:

``` r
devtools::load_all(".")
browseVignettes("MMRFBiolinks")
```

Get the list of the example data sets

``` r
data(package = "MMRFBiolinks")
```
