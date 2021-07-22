#source: drugCom
#install dependencies
#install.packages("devtools")
#devtools::install_github(c("RDocTaskForce/testextra","RDocTaskForce/parsetools","halpo/purrrogress"))
#install.packages(c("BiocManager","shiny","markdown","knitr","sweave","xtable","DT","httpuv","sourcetools"))
#BiocManager::install(c("EDASeq","genefilter","sva","limma","GenomicFeatures","GenomeInfoDb","GenomicRanges","SummarizedExperiment","EnsDb.Hsapiens.v79","S4Vectors","biomaRt","BiocStyle","edgeR","IRanges","TCGAbiolinks"))
#devtools::install_github("marziasettino/MMRFBiolinks")
#download and process
library(shiny)
#library(SummarizedExperiment)
#library(dplyr)
#library(DT)
#library(ggplot2)
#library(TCGAbiolinks)
#library(MMRFBiolinks)
#MMRFclin <- MMRFGDC_QueryClinic(type = "clinical")
#listSamples <- MMRFclin$bcr_patient_barcode
#query <- GDCquery(project = "MMRF-COMMPASS",data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",experimental.strategy = "RNA-Seq",workflow.type="HTSeq - FPKM",barcode = listSamples)
#GDCdownload(query, method = "api", files.per.chunk = 10)
#MMRnaseqSE <- MMRFGDC_prepare(query,save = TRUE ,save.filename = "GDCdata.rda",directory = "GDCdata",summarizedExperiment = TRUE)
#MMRFdataPrepro <- TCGAanalyze_Preprocessing(MMRnaseqSE)
#save(file="MMRFdataPrepro.rds",MMRFdataPrepro)
#save(file="MMRFclin.rds",MMRFclin)
#savehistory("R.history")
#load("MMRFclin.rds")
clinMMGateway<-MMRFBiolinks::clinMMGateway
#download.file("https://studntnu-my.sharepoint.com/:u:/g/personal/animeshs_ntnu_no/EfQcUVcV4m5IlA7YY9h6zFsBJwZhomsv-0r5uIgnd8Wt0g?e=UZC0vP","MMRFdataPrepro.rds")
#load("plotKM/MMRFdataPrepro.rds")
ui <- fluidPage(titlePanel("Drugs studied for Subjects with expression+clinical data in MMRF-COMMPASS"),textInput("DRUG",value="Dexamethasone",label = "Drug Name:",placeholder = "Daratumumab"),mainPanel(plotOutput("distPlot"),textOutput("textPlot")))
#MMRFRG_TimeBorPlot(clinMMGateway,"Dexamethasone","days")
#TCGAbiolinks::MMRFRG_TimeBorPlot(clinMMGateway,"Daratumumab","days")
#drugN<-"Dexamethasone"
#drugN<-"Daratumumab"
#MMRFclinDrug<- MMRFclin[which(unlist(MMRFclin[,"treatments"]) %in% c(drugN)==TRUE),]
#MMRFclinDrug<- MMRFclinDrug[!is.na(MMRFclinDrug$bcr_patient_barcode),]
#TCGAbiolinks::TCGAanalyze_SurvivalKM(MMRFclin,MMRFclin,Genelist = c("ENSG00000196449"),Survresult = T,ThreshTop=0.76,ThreshDown=0.33) #	Q86U90 (YRDC_HUMAN)
