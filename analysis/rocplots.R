modelROC <- function(data, model, allcalls) {
  vals <- (1:99)/100
  test.prediction <- predict(model, newdata=data)
  allcalls.prediction <- predict(model, newdata=allcalls)
  
  test.df <- data.frame(concordance=data$concordance, sample=data$sample, model=as.vector(test.prediction), validate_true=data$validate_true)
  allcalls.df <- data.frame(concordance=allcalls$concordance, sample=allcalls$sample, model=as.vector(allcalls.prediction))

  results <- lapply(vals, function(x) {
                            test.df$model <- as.vector(test.prediction > x);
                            allcalls.df$model <- as.vector(allcalls.prediction > x);
                            df <- corrected.accuracies(test.df, allcalls.df, "model");
                            df$thresh <- x;
                            return(df)})
  
  do.call(rbind,results)
}

library(randomForest)
rf <- function(formula, data) {
  newformula <- gsub('validate_true','factor(validate_true)', formula)
  return(randomForest(as.formula(newformula), data=data))
}

rfROC <- function(data, model) {
  prediction <- predict(model, data, type="prob")
  prediction <- as.vector(prediction[,2])
  vals <- (10:95)/100
  results <- sapply(vals, function(x) accuracies(data$validate_true, prediction > x))
  df <- data.frame(t(matrix(unlist(results),ncol=length(vals))))
  colnames(df) <- c("sensitivity", "precision", "f1")
  df <- df[,1:3]
  df$thresh <- vals
  return(df)
}

const.f1.rocs <- function(f1s) {
  sensrange <- function(f1) { mins <- max(f1/2+.01, f1/(2.-f1)+.01); seq(mins, .99, by=.01)}
  precrange <- function(sens, f1) {sens*f1/(2.*sens-f1)}
  sensprecdf <- function(f1) { s <- sensrange(f1); 
                               n <- length(s); 
                               data.frame(f1=rep(f1,n), sensitivity=s, precision=precrange(s, f1))}
  
  do.call(rbind, lapply(f1s, sensprecdf))
}

library(party)
library(ggplot2)
rocplot <- function(data, allcalls, formulae, names, callers, derived, title) {
  s <- corrected.accuracies.by.caller(data, allcalls, callers)
  s$derived <- derived
  
  allmodels <- data.frame()
  for (i in 1:length(names)) {
    for (trial in 1:10) {
      l <- split(data, sample(1:2, nrow(data), replace=TRUE))
      test <- l[[1]]; train <- l[[2]]
      
      treemodel <- ctree(as.formula(formulae[i]), data=train)
      results <- modelROC(test, treemodel, allcalls)
      results$model <- names[i]
      results$type <- "decision tree"
      
      allmodels <- rbind(allmodels, results)
      
      glmmodel <- glm(as.formula(formulae[i]), data=train)
      results <- modelROC(test, glmmodel, allcalls)
      results$model <- names[i]
      results$type <- "logistic regression"
      
      allmodels <- rbind(allmodels, results)
      
      #rfmodel <- rf(formulae[i], train)
      #results <- rfROC(test, rfmodel)
      #results$model <- names[i]
      #results$type <- "random forest"
    #  
      #allmodels <- rbind(allmodels, results)
    }
  }
  
  f1s <- seq(.4,.95,.05)
  roc.contours <- const.f1.rocs(f1s)
  ggplot() + geom_point(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller)) +
    geom_point(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller), color='blue') +
    geom_text(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller), angle=45, hjust=1, vjust=1) + 
    geom_text(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller, angle=+45), hjust=1, vjust=1, color='blue') + 
    geom_smooth(data=allmodels, aes(x=sensitivity, y=precision, color=type, linetype=model)) +
    geom_line(data=roc.contours, aes(x=sensitivity, y=precision, group=f1), color='grey', alpha=0.3) +
    ggtitle(title)
  
}
