modelROC <- function(data, model, allcalls, corrected=TRUE) {
  vals <- (5:45)/50
  test.prediction <- predict(model, newdata=data)
  allcalls.prediction <- predict(model, newdata=allcalls)
  
  test.df <- data.frame(concordance=data$concordance, sample=data$sample, model=as.vector(test.prediction), validate_true=data$validate_true)
  allcalls.df <- data.frame(concordance=allcalls$concordance, sample=allcalls$sample, model=as.vector(allcalls.prediction))

  accuracy.fun <- corrected.accuracies
  if (!corrected)
    accuracy.fun <- uncorrected.accuracies
  results <- lapply(vals, function(x) {
                            test.df$model <- as.vector(test.prediction > x);
                            allcalls.df$model <- as.vector(allcalls.prediction > x);
                            df <- accuracy.fun(test.df, allcalls.df, "model");
                            df$thresh <- x;
                            return(df)})
  
  do.call(rbind,results)
}

library(party)
library(ggplot2)
rocplot <- function(data, allcalls, formulae, names, callers, derived, title, 
                    xlim=c(0,1), ylim=c(0,1), include.models=TRUE, ntrials=10) {
  # limit samples in allcalls to those in data  
  data$sample <- factor(data$sample)
  allcalls <- allcalls[allcalls$sample %in% levels(data$sample), ]
  allcalls$sample <- factor(allcalls$sample)
  
  # get precision/recall values for callers
  s <- corrected.accuracies.by.caller(data, allcalls, callers)
  s$derived <- derived

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
  f1.df <- data.frame(sensitivity=sensitivity, precision=precision, f1=2*sensitivity*precision/(sensitivity+precision+1.e-7))
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

test.unseen.samples <- function(data, allcalls, formula, name, callers, derived, title, 
                          xlim=c(0,1), ylim=c(0,1)) {
  
  samples <- levels(data$sample)
  allcalls <- allcalls[allcalls$sample %in% samples, ]
  
  allmodels <- data.frame()
  
  ntrials <- as.integer(length(samples)/3)
  for (i in 1:ntrials) {
    samples.excluded <- samples[((i-1)*3+1):(i*3)]
    l <- split(data, data$sample %in% samples.excluded)
    seendata <- l[[1]]
    unseendata <- l[[2]]
    
    l <- split(seendata, sample(1:3, nrow(data), replace=TRUE)) 
    train <- rbind(l[[1]],l[[3]])
    test <- l[[2]][sample(1:nrow(l[[2]]), size=nrow(unseendata), replace=TRUE), ]
    
    treemodel <- ctree(as.formula(formula), data=train)
    results <- modelROC(test, treemodel, allcalls, FALSE)
    results$model <- name
    results$type <- "Decision tree"
    results$data <- "Seen samples"
    results$trial <- i
    allmodels <- rbind(allmodels, results)
    
    results <- modelROC(unseendata, treemodel, allcalls, FALSE)
    results$model <- name
    results$type <- "Decision tree"
    results$data <- "Unseen samples"
    results$trial <- i
    
    allmodels <- rbind(allmodels, results)
    
    glmmodel <- glm(as.formula(formula), data=train)
    results <- modelROC(test, glmmodel, allcalls, FALSE)
    results$model <- name
    results$type <- "Logistic regression"
    results$data <- "Seen samples"
    results$trial <- i
    
    allmodels <- rbind(allmodels, results)
    
    results <- modelROC(unseendata, glmmodel, allcalls, FALSE)
    results$model <- name
    results$type <- "Logistic regression"
    results$data <- "Unseen samples"
    results$trial <- i
    allmodels <- rbind(allmodels, results)
    
  }
  allmodels$data <- factor(allmodels$data)
  allmodels$model <- factor(allmodels$model)
  
  ggplot(allmodels) + geom_smooth(aes(x=sensitivity,y=precision,color=type,linetype=data))
  
  return(allmodels)
}
  
corrected.accuracies.by.caller.by.vaf.with.models <- function(data, allcalls, callers, derived, 
                                                       formula, name) {
  accuracies <- corrected.accuracies.by.caller.by.vaf(data, allcalls, callers)
  accuracies$derived <- ifelse(accuracies$caller %in% callers[derived], TRUE, FALSE)
  
  allcalls.df <- data.frame(concordance=allcalls$concordance, sample=allcalls$sample, binned_wgs_tvaf=allcalls$binned_wgs_tvaf)
  
  # make a train subdataset
  l <- split(data, sample(1:2, nrow(data), replace=TRUE))
  train <- l[[1]]; test <- l[[2]]
  
  treemodel <- ctree(as.formula(formula), data=train)
  glmmodel <- glm(as.formula(formula), data=train)
  
  thresh <- 0.5
  
  test.prediction <- predict(treemodel, newdata=test) > thresh
  test.df <- data.frame(concordance=test$concordance, sample=test$sample, decision_tree=as.vector(test.prediction), validate_true=test$validate_true, binned_wgs_tvaf=test$binned_wgs_tvaf)
  allcalls.df$decision_tree <- predict(treemodel, newdata=allcalls) > thresh
  
  treeaccuracies <- corrected.accuracies.by.caller.by.vaf(test.df, allcalls.df, c("decision_tree"))
  treeaccuracies$derived <- TRUE
  
  test.prediction <- predict(glmmodel, newdata=test) > thresh
  test.df$logistic_regression <- as.vector(test.prediction)
  allcalls.df$logistic_regression <- predict(treemodel, newdata=allcalls) > thresh
  
  glmaccuracies <- corrected.accuracies.by.caller.by.vaf(test.df, allcalls.df, c("logistic_regression"))
  glmaccuracies$derived <- TRUE
  
  rbind(accuracies,treeaccuracies,glmaccuracies)
}
                    
  
roc.plot.pdfs <- function() {
  
  rocplot(snvs, snv_calls, formulae, names, snv_callers, rep(FALSE, length(snv_callers)), 
          "SNV Calls: Array 2+3+4: Bin-Corrected", 
          xlim=c(.5,1), ylim=c(.7,1), include.models=FALSE)
  ggsave("plots/results/snv_roc_callers_only.pdf", width=10, height=10)

  rocplot(snvs, snv_calls, formulae, names, snv_callers_plus_derived, snv_derived, 
          "SNV Calls: Array 2+3+4: Bin-Corrected", 
          xlim=c(.5,1), ylim=c(.7,1), include.models=FALSE)
  ggsave("plots/results/snv_roc_callers_plus_derived.pdf", width=10, height=10)  
  
  rocplot(snvs, snv_calls, formulae, names, snv_callers_plus_derived, snv_derived, 
          "SNV Calls: Array 2+3+4: Bin-Corrected", 
          xlim=c(.5,1), ylim=c(.7,1), include.models=TRUE)
  ggsave("plots/results/snv_roc_callers_plus_derived_plus_models.pdf", width=12.7, height=10)  
  
  bad_indel_samples <- c("a34f1dba-5758-45c8-b825-1c888d6c4c13")
  goodindels <- indels[!indels$sample %in% bad_indel_samples, ]
  goodindel_calls <- indel_calls[!indel_calls$sample %in% bad_indel_samples, ]
  rocplot(indels, indel_calls, formulae, names, indel_callers_plus_derived, indel_derived, 
          "Indel Calls: Array 2+3+4 - 1 Sample (a34f): Bin-Corrected", 
          include.models=TRUE)
  ggsave("plots/results/indel_roc_callers_plus_derived_plus_models.pdf", width=12.7, height=10)  
  
  rocplot(indels[indels$repeat_count < 5,], indel_calls[indel_calls$repeat_count<5, ], formulae, names, indel_callers_plus_derived, indel_derived, 
          "Indel Calls< No Repeats: Array 2+3+4: Bin-Corrected", 
          include.models=TRUE)
  ggsave("plots/results/indel_roc_no_repeats_callers_plus_derived_plus_models.pdf", width=12.7, height=10)  
   
  formulae <- c("validate_true ~ sanger + embl_delly + broad_mutect")
  names <- ("Core callers")
  
  bad_sv_samples <- c("5e4bbb6b-66b2-4787-b8ce-70d17bc80ba8","bf95e410-b371-406c-a192-391d2fce94b2", "0e90fb64-00b2-4b53-bbc7-df8182b84060")
  goodsvs <- svs[!svs$sample %in% bad_sv_samples, ]
  goodsv_calls <- sv_calls[!sv_calls$sample %in% bad_sv_samples, ]
  rocplot(goodsvs, goodsv_calls, formulae, names, sv_callers_plus_derived, sv_derived, 
          "SV Calls: Array 2+3+4 - 3 samples (5e4b,bf95,0e90): Bin-Corrected", 
          include.models=FALSE)
  ggsave("plots/results/sv_roc_callers_plus_derived.pdf", width=10, height=10)  
}

merged.by.vaf.plots <- function() {
  core_callers <- c("broad_mutect","dkfz","sanger")
  results <- corrected.accuracies.by.caller.by.vaf.with.models(snvs, snv_calls, snv_callers_plus_derived, 
                                                               snv_derived, core_callers_formula,
                                                               "Core Callers + VAFs")
  results <- results[results$caller %in% core_callers | results$derived, ]
  print(plot.corrected.sensitivity.by.vaf(results) + ggtitle("SNV Sensitivity vs VAF, Core + Merged"))
  ggsave("plots/results/snv_sensitivity_vs_vaf_core_plus_merged.pdf", width=16, height=10)
  
  print(plot.corrected.precision.by.vaf(results) + ggtitle("SNV Precision vs VAF, Core + Merged"))
  ggsave("plots/results/snv_precision_vs_vaf_core_plus_merged.pdf", width=16, height=10)
  
  results <- corrected.accuracies.by.caller.by.vaf.with.models(indels, indel_calls, indel_callers_plus_derived, 
                                                               indel_derived, core_callers_formula,
                                                               "Core Callers + VAFs")
  results <- results[results$caller %in% core_callers | results$derived, ]
  print(plot.corrected.sensitivity.by.vaf(results) + ggtitle("Indel Sensitivity vs VAF, Core + Merged"))
  ggsave("plots/results/indel_sensitivity_vs_vaf_core_plus_merged.pdf", width=16, height=10)
  
  print(plot.corrected.precision.by.vaf(results) + ggtitle("Indel Precision vs VAF, Core + Merged"))
  ggsave("plots/results/indel_precision_vs_vaf_core_plus_merged.pdf", width=16, height=10)
}

roc.plot.uniform.pdfs <- function() {
  
  bad_snv_samples <- c("ee770885-b07c-4237-ae57-6eb52111446d","97449717-88cf-4caf-b4f3-d70f1bf7097d","a34f1dba-5758-45c8-b825-1c888d6c4c13")
  goodsnvs <- snvs[!snvs$sample %in% bad_snv_samples, ]
  goodsnv_calls <- snv_calls[!snvs$sample %in% bad_snv_samples, ]
  rocplot(snvs, snv_calls, formulae, names, snv_callers, rep(FALSE, length(snv_callers)), 
          "SNV Calls: All - 3 Samples(9744,a34f,ee77): Bin-Corrected", 
          include.models=FALSE)
  ggsave("~/Desktop/snvs.pdf", width=10, height=10)
  
  bad_indel_samples <- c("a34f1dba-5758-45c8-b825-1c888d6c4c13")
  goodindels <- indels[!indels$sample %in% bad_indel_samples, ]
  goodindel_calls <- indel_calls[!indel_calls$sample %in% bad_indel_samples, ]
  rocplot(indels, indel_calls, formulae, names, indel_callers, rep(FALSE, length(indel_callers)), 
          "Indel Calls: All - 1 Sample (a34f): Bin-Corrected", include.models=FALSE)
  ggsave("~/Desktop/indels.pdf", width=10, height=10)  
  
  rocplot(indels[indels$repeat_count < 5,], indel_calls[indel_calls$repeat_count<5, ], formulae, names, indel_callers, rep(FALSE, length(indel_callers)), 
          "Indel Calls, Repeat count < 5: All: Bin-Corrected", include.models=FALSE)
  ggsave("~/Desktop/indels-norepeats.pdf", width=10, height=10)  
  
  bad_sv_samples <- c("5e4bbb6b-66b2-4787-b8ce-70d17bc80ba8","bf95e410-b371-406c-a192-391d2fce94b2", "0e90fb64-00b2-4b53-bbc7-df8182b84060","a34f1dba-5758-45c8-b825-1c888d6c4c13")
  goodsvs <- svs[!svs$sample %in% bad_sv_samples, ]
  goodsv_calls <- sv_calls[!sv_calls$sample %in% bad_sv_samples, ]
  rocplot(goodsvs, goodsv_calls, formulae, names, sv_callers, rep(FALSE, length(sv_callers)), 
          "SV Calls: All - 3 samples (5e4b,a34f,bf95,0e90): Bin-Corrected", 
          include.models=FALSE)
  ggsave("~/Desktop/svs.pdf", width=10, height=10)  
}