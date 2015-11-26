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

library(party)
library(ggplot2)
rocplot <- function(data, allcalls, formulae, names, callers, derived, title, 
                    xlim=c(0,1), ylim=c(0,1), include.models=TRUE, ntrials=5) {
  # get precision/recall values for callers
  s <- corrected.accuracies.by.caller(data, allcalls, callers)
  s$derived <- derived

  # limit samples in allcalls to those in data  
  data$sample <- factor(data$sample)
  allcalls <- allcalls[allcalls$sample %in% levels(data$sample), ]
  allcalls$sample <- factor(allcalls$sample)
  
  # Generate 
  if (include.models) {
    allmodels <- data.frame()
    for (i in 1:length(names)) {
      for (trial in 1:ntrials) {
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
        
      }
    }
    allmodels <- aggregate(allmodels, by=list(allmodels$thresh, allmodels$type, allmodels$model), FUN=mean, na.rm=TRUE)[,1:6]
    colnames(allmodels) <- c("thresh","type","model","sensitivity","precision","f1")
  }
  
  df1 <- 0.05
  sensitivity <- rep(seq(0,1,df1), 1/df1+1)
  precision <- as.vector(t(matrix(sensitivity,ncol=as.integer(1/df1+1))))
  f1.df <- data.frame(sensitivity=sensitivity, precision=precision, f1=2*sensitivity*precision/(sensitivity+precision))
  f1.annotations <- seq(.1,.8,.1)
  f1.labels <- sapply(f1.annotations, function(x) paste("F1 = ",x))
  f1.annotations.df <- data.frame(sensitivity=f1.annotations, precision=f1.annotations, text=f1.labels)
  
  p <- ggplot() + geom_point(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller)) +
    geom_text(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller), angle=45, hjust=1, vjust=1) + 
    geom_contour(data=f1.df, aes(x=sensitivity, y=precision, z=f1), color='grey', alpha=0.5, binwidth=0.1) +
    geom_text(data=f1.annotations.df, aes(x=sensitivity, y=precision, label=text), color='grey', alpha=0.5) +
    theme(text = element_text(size=20)) +
    xlim(xlim) + ylim(ylim) + xlab('recall') +
    ggtitle(title)
  
  if (any(derived)) {
    p <- p + geom_point(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller), color='blue') +
      geom_text(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller, angle=+45), hjust=1, vjust=1, color='blue')
  }
  if (include.models) 
    p <- p + geom_line(data=allmodels, aes(x=sensitivity, y=precision, color=type, linetype=model))
                       
  p
}

roc.plot.pdfs <- function() {
  formulae <- c("validate_true ~ sanger + dkfz + broad_mutect + wgs_tvaf + wgs_nvaf + repeat_count")
  names <- ("Core callers + VAFs")
  
  rocplot(snvs, snv_calls, formulae, names, snv_callers, rep(FALSE, length(snv_callers)), 
          "SNV Calls: Array 2+3+4: Corrected Accuracies", 
          xlim=c(.5,1), ylim=c(.7,1), include.models=FALSE)
  ggsave("plots/results/snv_roc_callers_only.pdf", width=10, height=10)

  rocplot(snvs, snv_calls, formulae, names, snv_callers_plus_derived, snv_derived, 
          "SNV Calls: Array 2+3+4: Corrected Accuracies", 
          xlim=c(.5,1), ylim=c(.7,1), include.models=FALSE)
  ggsave("plots/results/snv_roc_callers_plus_derived.pdf", width=10, height=10)  
  
  rocplot(snvs, snv_calls, formulae, names, snv_callers_plus_derived, snv_derived, 
          "SNV Calls: Array 2+3+4: Corrected Accuracies", 
          xlim=c(.5,1), ylim=c(.7,1), include.models=TRUE)
  ggsave("plots/results/snv_roc_callers_plus_derived_plus_models.pdf", width=12.7, height=10)  
  
  rocplot(indels, indel_calls, formulae, names, indel_callers_plus_derived, indel_derived, 
          "Indel Calls: Array 2+3+4: Corrected Accuracies", 
          include.models=TRUE)
  ggsave("plots/results/indel_roc_callers_plus_derived_plus_models.pdf", width=12.7, height=10)  
}