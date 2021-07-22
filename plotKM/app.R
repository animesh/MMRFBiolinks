#source: plotKM
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
load("MMRFclin.rds")
#download.file("https://studntnu-my.sharepoint.com/:u:/g/personal/animeshs_ntnu_no/EfQcUVcV4m5IlA7YY9h6zFsBJwZhomsv-0r5uIgnd8Wt0g?e=UZC0vP","MMRFdataPrepro.rds")
load("MMRFdataPrepro.rds")
#legend("bottomleft")
#edit(TCGAbiolinks::TCGAanalyze_SurvivalKM,file="bottomleft.r")
#source("bottomleft.r")
#body(blKM)<-source("bottomleft.r")
#https://stackoverflow.com/a/26089312/1137129
blKM<-TCGAbiolinks::TCGAanalyze_SurvivalKM#<-TCGAanalyze_SurvivalKM
body(blKM)[[25]][[4]]<-substitute({
  cat(paste0((ngenes - i), "."))
  mRNAselected <- as.matrix(rownames(dataNormal))[i]
  mRNAselected_values <- dataCancer[rownames(dataCancer) == 
                                      mRNAselected, ]
  mRNAselected_values_normal <- dataNormal[rownames(dataNormal) == 
                                             mRNAselected, ]
  if (all(mRNAselected_values == 0)) 
    next
  tabSurv_Matrix[i, "mRNA"] <- mRNAselected
  mRNAselected_values_ordered <- sort(mRNAselected_values, 
                                      decreasing = TRUE)
  mRNAselected_values_ordered_top <- as.numeric(quantile(as.numeric(mRNAselected_values_ordered), 
                                                         ThreshTop)[1])
  mRNAselected_values_ordered_down <- as.numeric(quantile(as.numeric(mRNAselected_values_ordered), 
                                                          ThreshDown)[1])
  mRNAselected_values_newvector <- mRNAselected_values
  if (!is.na(mRNAselected_values_ordered_top)) {
    numberOfSamples <- length(mRNAselected_values_ordered)
    lastelementTOP <- max(which(mRNAselected_values_ordered > 
                                  mRNAselected_values_ordered_top))
    firstelementDOWN <- min(which(mRNAselected_values_ordered <= 
                                    mRNAselected_values_ordered_down))
    samples_top_mRNA_selected <- names(mRNAselected_values_ordered[1:lastelementTOP])
    samples_down_mRNA_selected <- names(mRNAselected_values_ordered[firstelementDOWN:numberOfSamples])
    samples_UNCHANGED_mRNA_selected <- names(mRNAselected_values_newvector[which((mRNAselected_values_newvector) > 
                                                                                   mRNAselected_values_ordered_down & mRNAselected_values_newvector < 
                                                                                   mRNAselected_values_ordered_top)])
    cfu_onlyTOP <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% 
                                  samples_top_mRNA_selected, ]
    cfu_onlyDOWN <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% 
                                   samples_down_mRNA_selected, ]
    cfu_onlyUNCHANGED <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% 
                                        samples_UNCHANGED_mRNA_selected, ]
    cfu_ordered <- NULL
    cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
    cfu <- cfu_ordered
    ttime <- as.numeric(cfu[, "days_to_death"])
    sum(status <- ttime > 0)
    deads_complete <- sum(status <- ttime > 0)
    ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
    deads_top <- sum(ttime_only_top > 0)
    if (dim(cfu_onlyDOWN)[1] >= 1) {
      ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
      deads_down <- sum(ttime_only_down > 0)
    }
    else {
      deads_down <- 0
    }
    tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
    tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
    tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
    tabSurv_Matrix[i, "Mean Normal"] <- mean(as.numeric(mRNAselected_values_normal))
    dataCancer_onlyTop_sample <- dataCancer[, samples_top_mRNA_selected, 
                                            drop = FALSE]
    dataCancer_onlyTop_sample_mRNASelected <- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == 
                                                                          mRNAselected, ]
    dataCancer_onlyDown_sample <- dataCancer[, samples_down_mRNA_selected, 
                                             drop = FALSE]
    dataCancer_onlyDown_sample_mRNASelected <- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == 
                                                                            mRNAselected, ]
    tabSurv_Matrix[i, "Mean Tumor Top"] <- mean(as.numeric(dataCancer_onlyTop_sample_mRNASelected))
    tabSurv_Matrix[i, "Mean Tumor Down"] <- mean(as.numeric(dataCancer_onlyDown_sample_mRNASelected))
    ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
    ttime[which(ttime == -Inf)] <- 0
    ttime <- survival::Surv(ttime, status)
    rownames(ttime) <- rownames(cfu)
    legendHigh <- paste(mRNAselected, "High")
    legendLow <- paste(mRNAselected, "Low")
    tabSurv_pvalue <- tryCatch({
      tabSurv <- survival::survdiff(ttime ~ c(rep("top", 
                                                  nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN))))
      tabSurv_chis <- unlist(tabSurv)$chisq
      tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), 
                                              df = 1))
    }, error = function(e) {
      return(Inf)
    })
    tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue
    if (Survresult == TRUE) {
      titlePlot <- paste("Kaplan-Meier Survival analysis, pvalue=", 
                         tabSurv_pvalue)
      plot(survival::survfit(ttime ~ c(rep("low", 
                                           nrow(cfu_onlyTOP)), rep("high", nrow(cfu_onlyDOWN)))), 
           col = c("green", "red"), main = titlePlot, 
           xlab = "Days", ylab = "Survival")
      legend("bottomleft",legend = c(legendLow, legendHigh), 
             col = c("green", "red"), text.col = c("green","red"), 
             pch = 15)
      print(tabSurv)
    }
  }
})
#https://stackoverflow.com/questions/19226816/how-can-i-view-the-source-code-for-a-function
#newDef <- deparse(TCGAbiolinks::TCGAanalyze_SurvivalKM)
#newDef[grep("                legend",newDef)]<- "                legend(100, 1, legend = c(legendLow, legendHigh), "
#blKM <- eval(parse(text=newDef))
ui <- fluidPage(titlePanel("Plot Kaplan–Meier estimator for Subjects with expression+clinical data in MMRF-COMMPASS"),textInput("ENSG",value="ENSG00000196976",label = "Ensemble Gene ID:",placeholder = "ENSG00000196976"),mainPanel(plotOutput("distPlot"),textOutput("textPlot")))
server <- function(input, output) {
#  output$distPlot <- renderPlot({TCGAbiolinks::TCGAanalyze_SurvivalKM(MMRFclin,MMRFdataPrepro,Genelist = input$ENSG,Survresult = T,ThreshTop=0.76,ThreshDown=0.33)})
#  output$textPlot <- renderPrint({TCGAbiolinks::TCGAanalyze_SurvivalKM(MMRFclin,MMRFdataPrepro,Genelist = input$ENSG,Survresult = F,ThreshTop=0.76,ThreshDown=0.33)})
  output$distPlot <- renderPlot({blKM(MMRFclin,MMRFdataPrepro,Genelist = input$ENSG,Survresult = T,ThreshTop=0.76,ThreshDown=0.33)})
  output$textPlot <- renderPrint({blKM(MMRFclin,MMRFdataPrepro,Genelist = input$ENSG,Survresult = F,ThreshTop=0.76,ThreshDown=0.33)})
}
shinyApp(ui = ui, server = server)
