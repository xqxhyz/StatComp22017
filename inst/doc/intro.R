## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
#load('../data/data.rda') #data
#load('../data/data2.rda') #data2
#cutpoint <- 0.0002820585
#SAIDP(data=data, data2=data2, cutpoint=0.0002820585)

## -----------------------------------------------------------------------------
library(e1071)

## -----------------------------------------------------------------------------
# fit_continuous(
#   data,
#   frac = 0.8,
#   models = c("LM","SVM"),
# )
fit_continuous <- function(data, frac=0.8, models=c("LM","SVM")){
  set.seed(516)
  # split train data and test data
  y = 0 
  m <- nrow(data)
  n <- ncol(data)
  colnames(data) <- c(paste0('v',1:(n-1)),'y')
  train_index <- sample(c(1:m), floor(frac*m))
  train_data <- data[train_index,]
  test_data <- data[-train_index,]
  
  MSE <- function(pred,true){mean((pred-true)^2)}
  MSEs <- c()
  row_names <- c()
  model = list()
  
  if("LM" %in% models){
    model_lm <- lm(y~., data = train_data)
    pred_lm <- predict(model_lm, test_data)
    MSE_lm <- MSE(pred_lm,test_data$y)
    plot(test_data$y, pred_lm, main="LM", xlab="Ground Truth", ylab="Prediction")
    abline(0, 1, col="orange")
    MSEs <- c(MSEs, MSE_lm)
    row_names <- c(row_names, "LM")
    model1 <- list(model_lm)
    model <- append(model, model1)
  }
  
  if("SVM" %in% models){
    model_svm <- svm(y~., train_data, type="eps-regression", kernel="linear")
    pred_svm <- predict(model_svm, test_data)
    MSE_svm <- MSE(pred_svm, test_data$y)
    plot(pred_svm, test_data$y, main="SVR", xlab="Ground Truth", ylab="Prediction")
    abline(0, 1, col="orange")
    MSEs <- c(MSEs, MSE_svm)
    row_names <- c(row_names, "SVM")
    model2 <- list(model_svm)
    model <- append(model, model2)
  }
  
  MSES <- data.frame(model=row_names, MSE=MSEs)
  return(list(MSES,model))
}

## -----------------------------------------------------------------------------
set.seed(516)
x1 <- 1:100
x2 <- seq(2, 18, length.out=100) + rnorm(100, 3, 4)
y <- 4*x1 - 3*x2 + rnorm(100, 0, 2)
data <- data.frame(x1=x1, x2=x2, y=y)
results <- fit_continuous(data, frac=0.8, c("LM","SVM"))
results[[1]]
results[[2]][1]

