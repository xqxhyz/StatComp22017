#' @title Train models for continuous label
#' @description This function can automatically split the train test datasets, and train LM and SVM models for continuous label and gain the MSE on test set and models. Also plot the prediction and ground truth of test data. 
#' @param data A data frame with feature and label.
#' @param frac The proportion of train_data. The default value is 0.8.
#' @param models A vector of model names want to be trained. The default value is c("LM","SVM").
#' @return A list of MSE and models.
#' @import e1071
#' @examples 
#' \dontrun{
#' set.seed(516)
#' x1 <- 1:100
#' x2 <- seq(2, 18, length.out=100) + rnorm(100, 3, 4)
#' y <- 4x1 - 3x2 + rnorm(100, 0, 2)
#' data <- data.frame(x1=x1, x2=x2, y=y)
#' results <- fit_continuous(data, frac=0.8, models=c("LM","SVM"))
#' }
#' @export
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