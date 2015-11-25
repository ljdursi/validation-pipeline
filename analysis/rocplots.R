modelROC <- function(data, model, allcalls) {
  vals <- (2:23)/25
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

library(party)
library(ggplot2)
rocplot <- function(data, allcalls, formulae, names, callers, derived, title, xlim=c(0,1), ylim=c(0,1)) {
  s <- corrected.accuracies.by.caller(data, allcalls, callers)
  s$derived <- derived
  
  data$sample <- factor(data$sample)
  allcalls <- allcalls[allcalls$sample %in% levels(data$sample), ]
  allcalls$sample <- factor(allcalls$sample)
  
  weights <- xtabs(~concordance+sample,data=data)/xtabs(~concordance+sample,data=allcalls)
  data$weights <- as.vector(apply(data[,c("concordance","sample")], 1, function(x) weights[as.integer(x[1]), x[2]]))
  data$weights[data$weights < 1] <- 1
  data$weights[data$weights > 10] <- 10
  data$weights <- as.integer(data$weights)
  
  allmodels <- data.frame()
  for (i in 1:length(names)) {
    for (trial in 1:10) {
      l <- split(data, sample(1:2, nrow(data), replace=TRUE))
      test <- l[[1]]; train <- l[[2]]
      
      treemodel <- ctree(as.formula(formulae[i]), data=train, weights=train$weights)
      results <- modelROC(test, treemodel, allcalls)
      results$model <- names[i]
      results$type <- "decision tree"
      
      allmodels <- rbind(allmodels, results)
      
      glmmodel <- glm(as.formula(formulae[i]), data=train, weights=train$weights)
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
  allmodels <- aggregate(allmodels, by=list(allmodels$thresh, allmodels$type, allmodels$model), FUN=mean, na.rm=TRUE)[,1:6]
  colnames(allmodels) <- c("thresh","type","model","sensitivity","precision","f1")
  
  df1 <- 0.05
  sensitivity <- rep(seq(0,1,df1), 1/df1+1)
  precision <- as.vector(t(matrix(sensitivity,ncol=as.integer(1/df1+1))))
  f1.df <- data.frame(sensitivity=sensitivity, precision=precision, f1=2*sensitivity*precision/(sensitivity+precision))
  f1.annotations <- seq(.1,.8,.1)
  f1.labels <- sapply(f1.annotations, function(x) paste("F1 = ",x))
  f1.annotations.df <- data.frame(sensitivity=f1.annotations, precision=f1.annotations, text=f1.labels)
  
  ggplot() + geom_point(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller)) +
    geom_point(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller), color='blue') +
    geom_text(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller), angle=45, hjust=1, vjust=1) + 
    geom_text(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller, angle=+45), hjust=1, vjust=1, color='blue') + 
    geom_line(data=allmodels, aes(x=sensitivity, y=precision, color=type, linetype=model)) +
    geom_contour(data=f1.df, aes(x=sensitivity, y=precision, z=f1), color='grey', alpha=0.5, binwidth=0.1) +
    geom_text(data=f1.annotations.df, aes(x=sensitivity, y=precision, label=text), color='grey', alpha=0.5) +
    theme(text = element_text(size=20)) +
    xlim(xlim) + ylim(ylim) + xlab('recall') +
    ggtitle(title)
}
