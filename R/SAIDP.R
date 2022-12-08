#' @title Survival analysis in digital pathology
#' @description This function can conduct data preprocessing for TIL related information obtained by the model and obtain a series of survival analysis results.
#' @param data Basic information of each patient in the dataset, such as survival time, survival state, age, etc.
#' @param data2 Information about TIL of each patient obtained by neural network model.
#' @param cutpoint Determine the size of the TIL Score based on the cutpoint.
#' @return A series of survival analyses of datasets.
#' @import DT
#' @import survival
#' @import survminer
#' @import ggplot2
#' @import DAAG
#' @import Rcpp
#' @import boot
#' @import bootstrap
#' @import knitr
#' @import latex2exp
#' @import microbenchmark
#' @import pander
#' @import rmarkdown
#' @import scales
#' @import stats
#' @import graphics
#' @examples 
#' \dontrun{
#' #load("data\\data.rda")
#' #load("data\\data2.rda")
#' #SAIDP(data=data, data2=data2, cutpoint=0.0002820585)
#' }
#' @export
SAIDP <- function(data, data2, cutpoint){
data$Survival.Time.months. = as.numeric(as.character(data$Survival.Time.months.))
data_surv <- Surv(data$Survival.Time.months., data$Vital.Status=="Dead")
kmfit1 <- survfit(data_surv~1)
ggsurvplot(kmfit1, data=data)

data2 <- cbind(data2,
               Survival.Time.months.=NA, Vital.Status=NA,
               Age=NA,
               TIL_status=NA)
data2 <- data2[!(is.na(data2$TIL)),]

for (i in 1:nrow(data2)){
  id <- data2$Patient_ID[i]
  imap <- which(data$Patient_ID==id)[1]
  data2$Survival.Time.months.[i] <- data$Survival.Time.months[imap]     
  if (data$Vital.Status[imap]=="Alive"){
    data2$Vital.Status[i] = 0
  }
  if (data$Vital.Status[imap]=="Dead"){
    data2$Vital.Status[i] = 1
  }
  data2$Age[i] <- data$Age[imap]
  data2$TIL_status = as.numeric((data2$TIL >= cutpoint))+1 #1???TIL??????2???TIL???
}

data2 <- data2[data2$TIL>0,]

res.cut_OS <- surv_cutpoint(data2, time = "Survival.Time.months.", event = "Vital.Status",
                            variables = c("TIL"))
summary(res.cut_OS)

a <- ggplot(data2, aes(x = TIL))
a + xlab("TIL Score")+ylab("Count")+geom_histogram(bins = 30, color = "black", fill = "gray" ) 

attach(data2)
data_surv_os <- Surv(Survival.Time.months., Vital.Status==1) #1:Dead;0:Alive

kmfit_os <- survfit(data_surv_os~TIL_status)

survdiff(data_surv_os~data2$TIL_status)
plot(kmfit_os, col=c(1,2), xlab="Time (months)", ylab="Survival probability", main="KM Curve")
legend("bottomleft", title = "TIL-score", c("high", "low"), col=c(1,2), lty=1)

ggsurvplot(kmfit_os, pval = TRUE, data = data2, title = "KM Curve", palette = "nature", legend.labs = c("Low Risk","High Risk"))+xlab("Time(month)")+ylab("survival probability")
}