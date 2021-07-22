library(shiny)
clinMMGateway<-MMRFBiolinks::clinMMGateway
server <- function(input, output) {
    output$distPlot <- renderPlot({barplot(summary(as.factor(clinMMGateway$trtname)),horiz = T)})
    output$textPlot <- renderPrint(summary(as.factor(clinMMGateway$trtname)))
    #MMRFRG_GetBorInfo(clinMMGateway)
    #MMRFRG_TimeBorPlot(clinMMGateway,"Dexamethasone","days")
}
#MMRFRG_TimeBorPlot(clinMMGateway,"Dexamethasone","days")
#MMRFRG_TimeBorPlot(clinMMGateway,"Daratumumab","days")
#MMRFclinDrug<- MMRFclin[which(unlist(MMRFclin[,"treatments"]) %in% c("Dexamethasone")==TRUE),]
#MMRFclinDrug<- MMRFclinDrug[!is.na(MMRFclinDrug$bcr_patient_barcode),]
#TCGAbiolinks::TCGAanalyze_SurvivalKM(MMRFclin,MMRFdataPrepro,Genelist = c("ENSG00000196449"),Survresult = T,ThreshTop=0.66,ThreshDown=0.33) #	Q86U90 (YRDC_HUMAN)
#TCGAbiolinks::TCGAanalyze_SurvivalKM(MMRFclin,MMRFclinDrug,Genelist = c("ENSG00000196449"),Survresult = T,ThreshTop=0.76,ThreshDown=0.33) #	Q86U90 (YRDC_HUMAN)
