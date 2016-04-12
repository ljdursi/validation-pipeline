dummy.lines <- function(template, n, vafprobs=c(1,0,0,0,0)) {
  
  tmp <- template[1:n,]
  tmp$chrom <- "chr0"
  tmp$val_tvaf <- -1
  tmp$val_nvaf <- -1
  tmp$val_tdepth <- -1
  tmp$val_ndepth <- -1
  tmp$wgs_tvaf <- -1
  tmp$wgs_nvaf <- -1
  tmp$binned_wgs_tvaf <- sample(levels(template$binned_wgs_tvaf), n, replace=TRUE, prob=vafprobs)

  tmp$repeat_count <- -1
  
  tmp$concordance <- 0
  tmp$adiscan <- 0
  tmp$broad_mutect <- 0
  tmp$dkfz <- 0
  tmp$lohcomplete <- 0
  tmp$mda_hgsc_gatk_muse <- 0
  tmp$oicr_bl <- 0
  tmp$oicr_sga <- 0
  tmp$sanger <- 0
  tmp$smufin <- 0
  tmp$wustl <- 0
  tmp$common_sample <- TRUE
  tmp$union <- 0
  tmp$intersect2 <- 0
  tmp$intersect3 <- 0
  
  tmp
}

made.up.validation.calls <- function(template, passall=TRUE, privateaccuracies=NA, npersample=50, callername="new_mutect", vafprobs=c(1,0,0,0,0)) {
  sample_levels <- levels(template$sample)
  result <- data.frame()
  
  for (isample in seq(length(levels(template$sample)))) {
    sample <- sample_levels[isample]
    tmp <- dummy.lines(template, npersample, vafprobs)
    
    tmp$sample <- sample
    tmp[[callername]] <- 1
    
    if (!is.na(privateaccuracies)) {
      frac <- privateaccuracies[sample,]$precision
      tmp$validate_true <- runif(npersample) > frac
      tmp$status <- ifelse(tmp$validate_true, "PASS", "GERMLINE")
    } else if (passall) {
      tmp$status <- "PASS"
      tmp$validate_true <- TRUE
    } else {
      tmp$status <- "NOTSEEN"
      tmp$validate_true <- FALSE
    }
    
    result <- rbind(result, tmp)
  }
  result
}

estimates <- function(snvs, snv_calls,
                      filename='~/Desktop/otherbroad/new_mutect.csv',
                      newcaller='new_mutect',
                      originalcaller='broad_mutect',
                      shortcaller='mutect') {
  
  privatecaller <- paste0(shortcaller, '_private')
  newcaller <- paste0('new_',shortcaller)
  
  snvs[[privatecaller]] <- snvs[[originalcaller]] & snvs$concordance == 1
  snv_calls[[privatecaller]] <- snv_calls[[originalcaller]] & snv_calls$concordance == 1
  
  old_mutect_private_acc_by_sample <- corrected.accuracies(snvs, snv_calls, privatecaller, combine=c('concordance'))
  results <- corrected.accuracies(snvs, snv_calls, originalcaller, combine=c('concordance'))
  
  results$caller <- paste0(shortcaller,"_original")
  results<- melt(results, id=c('sample','caller'))
  
  dkfz <- corrected.accuracies(snvs, snv_calls, 'dkfz', combine=c('concordance'))
  sanger <- corrected.accuracies(snvs, snv_calls, 'sanger', combine=c('concordance'))
  broad <- corrected.accuracies(snvs, snv_calls, 'broad_mutect', combine=c('concordance'))
  others <- rbind(dkfz, sanger, broad)
  others <- melt(others, id=c('sample','caller'))
  results <- rbind(results, others)
  
  new_mutect <- read.csv(filename)
  new_mutect <- new_mutect[new_mutect$sample %in% levels(snvs$sample),]
  new_mutect$ref <- as.character(new_mutect$ref)
  new_mutect$alt <- as.character(new_mutect$alt)
  
  # validation data; we're not going to change these, except add a 'new_mutect' caller to existing calls where appropriate
  new_snvs <- merge(snvs, new_mutect, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snvs[[newcaller]][is.na(new_snvs[[newcaller]])] <- 0

  # all calls - include the novel mutect calls, with concordance 0:
  new_snv_calls <- merge(snv_calls, new_mutect, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=TRUE)
  new_snv_calls[[newcaller]][is.na(new_snv_calls[[newcaller]])] <- 0
  new_snv_calls$concordance[is.na(new_snv_calls$concordance)] <- 0
  
  # accuracies assuming all novel calls are wrong
  new_snvs_with_novel_false <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=FALSE, callername=newcaller))
  accuracies <- corrected.accuracies(new_snvs_with_novel_false, new_snv_calls, newcaller, combine=c('concordance'))
  accuracies$caller <- paste0(shortcaller,"_novel_calls_bad")
  accuracies <- melt(accuracies, id=c('sample','caller'))
  results <- rbind(results, accuracies)
  
  # accuracies assuming all novel calls are right
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=TRUE, callername=newcaller))
  accuracies <- corrected.accuracies(new_snvs_with_novel_true, new_snv_calls, newcaller, combine=c('concordance'))
  accuracies$caller <- paste0(shortcaller,"_novel_calls_good")
  accuracies <- melt(accuracies, id=c('sample','caller'))
  results <- rbind(results, accuracies)
  
  # accuracies extrapolating the per-sample private call accuracy
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, callername=newcaller, privateaccuracies=old_mutect_private_acc_by_sample))
  accuracies <- corrected.accuracies(new_snvs_with_novel_true, new_snv_calls, newcaller, combine=c('concordance'))
  accuracies$caller <- paste0(shortcaller,"_novel_calls_interpolated")
  accuracies <- melt(accuracies, id=c('sample','caller'))
  results <- rbind(results, accuracies)

  # sort samples by number of intersect2 calls:
  sorted_samples <- names(sort(tapply(snv_calls$intersect2, snv_calls$sample, FUN=sum)))
  results$sample <- factor(results$sample, levels=sorted_samples)
  
  results$base_caller <- as.character(results$caller)
  results$base_caller[grep(shortcaller, results$base_caller)] <- shortcaller
  results$base_caller <- as.factor(results$base_caller)
  results$estimated <- FALSE
  results$estimated[results$base_caller == shortcaller] <- TRUE
  results$estimated[results$caller == paste0(shortcaller,"_original")] <- FALSE
  results
}

estimates.by.vaf <- function(snvs, snv_calls, snv_callers,
                             filename='~/Desktop/otherbroad/new_mutect.csv',
                             originalcaller='broad_mutect',
                             shortcaller='mutect',
                             vafprobs=c(1,0,0,0,0)) {
  
  privatecaller <- paste0(shortcaller, '_private')
  newcaller <- paste0('new_',shortcaller)
  
  snvs[[privatecaller]] <- snvs[[originalcaller]] & snvs$concordance == 1
  snv_calls[[privatecaller]] <- snv_calls[[originalcaller]] & snv_calls$concordance == 1
  old_caller_private_acc_by_sample <- corrected.accuracies(snvs, snv_calls, privatecaller, combine=c('concordance'))
  
  results <- corrected.accuracies.by.caller.by.vaf(snvs, snv_calls, snv_callers)

  new_calls <- read.csv(filename)
  new_calls <- new_calls[new_calls$sample %in% levels(snvs$sample),]
  new_calls$ref <- as.character(new_calls$ref)
  new_calls$alt <- as.character(new_calls$alt)
  
  # validation data; we're not going to change these, except add a 'new_mutect' caller to existing calls where appropriate
  new_snvs <- merge(snvs, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snvs[[newcaller]][is.na(new_snvs[[newcaller]])] <- 0
  
  # all calls - include the novel mutect calls, with concordance 0:
  new_snv_calls <- merge(snv_calls, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=TRUE)
  new_snv_calls[[newcaller]][is.na(new_snv_calls[[newcaller]])] <- 0
  new_snv_calls$concordance[is.na(new_snv_calls$concordance)] <- 0
  new_snv_calls$binned_wgs_tvaf[is.na(new_snv_calls$binned_wgs_tvaf)] <- sample(levels(snv_calls$binned_wgs_tvaf),
                                                                                sum(is.na(new_snv_calls$binned_wgs_tvaf)),
                                                                                replace=TRUE,
                                                                                prob=vafprobs)
  
  # accuracies assuming all novel calls are wrong
  new_snvs_with_novel_false <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=FALSE, callername=newcaller, vafprobs=vafprobs))
  accuracies <- corrected.accuracies.by.caller.by.vaf(new_snvs_with_novel_false, new_snv_calls, newcaller)
  accuracies$caller <- paste0(shortcaller,"_novel_calls_bad")
  results <- rbind(results, accuracies)
  
  # accuracies assuming all novel calls are right
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, passall=TRUE, callername=newcaller, vafprobs=vafprobs))
  accuracies <- corrected.accuracies.by.caller.by.vaf(new_snvs_with_novel_true, new_snv_calls, newcaller)
  accuracies$caller <- paste0(shortcaller, "_novel_calls_good")
  results <- rbind(results, accuracies)
  
  # accuracies extrapolating the per-sample private call accuracy
  new_snvs_with_novel_true <- rbind(new_snvs, made.up.validation.calls(new_snvs, callername=newcaller, privateaccuracies=old_caller_private_acc_by_sample, vafprobs=vafprobs))
  accuracies <- corrected.accuracies.by.caller.by.vaf(new_snvs_with_novel_true, new_snv_calls, newcaller)
  accuracies$caller <- paste0(shortcaller, "_novel_calls_iterpolated")
  results <- rbind(results, accuracies)
  
  results <- melt(results, id=c('VAF','caller'))
  results$base_caller <- as.character(results$caller)
  results$base_caller[grep(shortcaller, results$base_caller)] <- shortcaller
  results$estimated <- FALSE
  results$estimated[results$base_caller == shortcaller] <- TRUE
  results$estimated[results$caller == originalcaller] <- FALSE
  results
}

muse.estimates <- function(snvs, snv_calls, filename='~/Desktop/othermuse/new_muse.csv') {
  estimates(snvs, snv_calls, filename=filename, newcaller='new_muse',
            originalcaller='mda_hgsc_gatk_muse', shortcaller='muse')
}

muse.vaf.estimates <- function(snvs, snv_calls, snv_callers, 
                               filename='~/Desktop/newmuse/new_muse.csv', 
                               origcaller='mda_hgsc_gatk_muse', 
                               shortcaller='muse', 
                               vafprobs=c(0.5,0.5,0,0,0)) {
  estimates.by.vaf(snvs, snv_calls, snv_callers, filename=filename, originalcaller=origcaller, shortcaller=shortcaller, vafprobs=vafprobs)
}

add.newcaller.as.feature <- function(snvs, snv_calls, newcaller='new_muse', filename='~/Desktop/newmuse/new_muse.csv') {
  
  new_calls <- read.csv(filename)
  new_calls <- new_calls[new_calls$sample %in% levels(snvs$sample),]
  new_calls$ref <- as.character(new_calls$ref)
  new_calls$alt <- as.character(new_calls$alt)
  
  # validation data; we're not going to change these, except add a 'newcaller' caller to existing calls where appropriate
  new_snvs <- merge(snvs, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snvs[[newcaller]][is.na(new_snvs[[newcaller]])] <- 0
  
  # all calls - include the novel mutect calls, with concordance 0:
  new_snv_calls <- merge(snv_calls, new_calls, by=c("chrom","pos","alt","sample"), all.x=TRUE, all.y=FALSE)
  new_snv_calls[[newcaller]][is.na(new_snv_calls[[newcaller]])] <- 0
  
  return(list(new_snvs=new_snvs, new_snv_calls=new_snv_calls))
}

library(reshape2)
library(ggplot2)

boxplots <- function(snvs, snv_calls, callers=c('broad_mutect','dkfz','sanger'), samples=levels(snvs$sample)) {
  result.list <- lapply(callers, function(x) corrected.accuracies(snvs, snv_calls, x, combine=c('concordance')))
  results <- do.call(rbind, result.list)
  results <- results[results$sample %in% samples,]
  results <- melt(results,id=c('caller','sample'))
  results = results[results$value != 0,]
  ggplot(results) + geom_boxplot(aes(x=caller,y=value,color=caller)) + facet_grid(variable ~ .) + ggtitle("Distribution of Per-Sample Accuracies") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


applymodel <- function(snvs, snv_calls, learn_model, predict_model, formula, model_name, samples=levels(snvs$sample), seed=NA) {
  if (!is.na(seed)) set.seed(seed)

  copy_snvs <- snvs
  copy_snvs$row <- 1:nrow(snvs)

  copy_snv_calls <- snv_calls
  copy_snv_calls$row <- 1:nrow(snv_calls)

  copy_snvs <- copy_snvs[copy_snvs$union == 1,]
  copy_snv_calls <- copy_snv_calls[copy_snv_calls$union == 1,]

  nsamples <- length(samples)
  train_idxs <- sample(1:nsamples, nsamples/2, replace=FALSE)
  train_samples <- samples[train_idxs]
  test_samples <- samples[-train_idxs]
  
  train_snvs <- copy_snvs[copy_snvs$sample %in% train_samples,]
  train_snv_calls <- copy_snv_calls[copy_snv_calls$sample %in% train_samples, ]
  test_snvs <- copy_snvs[copy_snvs$sample %in% test_samples, ]
  test_snv_calls <- copy_snv_calls[copy_snv_calls$sample %in% test_samples, ]
  
  # learn on the train and apply to test
  model <- learn_model(formula, train_snvs)
  test_snv_calls_predictions <- predict_model(model, formula, test_snv_calls)
  test_snvs_predictions <- predict_model(model, formula, test_snvs)
  
  # and vice versa, so that all samples have equally good/bad predictions
  model <- learn_model(formula, test_snvs)
  train_snv_calls_predictions <- predict_model(model, formula, train_snv_calls)
  train_snvs_predictions <- predict_model(model, formula, train_snvs)
  
  copy_snvs <- snvs
  copy_snvs[[model_name]] <- 0
  copy_snvs[[model_name]][test_snvs$row] <- test_snvs_predictions
  copy_snvs[[model_name]][train_snvs$row] <- train_snvs_predictions
  
  copy_snv_calls <- snv_calls
  copy_snv_calls[[model_name]] <- 0
  copy_snv_calls[[model_name]][test_snv_calls$row] <- test_snv_calls_predictions
  copy_snv_calls[[model_name]][train_snv_calls$row] <- train_snv_calls_predictions
  
  return(list(snv_calls=copy_snv_calls, snvs=copy_snvs))
}

glmlearn <- function(formula, data) { glm(formula, data, family="binomial")}
glmpredict <- function(model, formula, data) { predict(model, data, type="response")}

ctreepredict <- function(model, formula, data) { predict(model, newdata=data) }

require(glmnet)
require(party)

glmnetlearn <- function(formula, data) {
  mydata <- data[complete.cases(data),]
  X <- model.matrix(formula, data=mydata)
  cv.glmnet(X, mydata$validate_true, family="binomial", foldid=as.integer(droplevels(mydata$sample)))
}

glmnetpredict <- function(model, formula, data) {
  results <- rep(0,nrow(data))
  X <- model.matrix(formula, data=data)
  used_rows <- as.integer(dimnames(X)[[1]])
  results[used_rows] <- predict(model, X, type="response", s="lambda.1se")[,1]
  return(results)
}

binarize <- function(data, thresh=0.5) {
    ifelse(data > thresh, 1, 0)
}


simple_indel_models <- function(validated.calls, all.calls, seed=1) {
  l <-  applymodel(validated.calls, all.calls, glmlearn, glmpredict, 
                   as.formula("validate_true ~ wgs_nvaf + wgs_tvaf + varlen + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger"), 
                   "logistic_regression", seed=1)
  
  l <-  applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
                   as.formula("validate_true ~ wgs_nvaf + wgs_tvaf + varlen+ gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger"), 
                   "decision_tree", seed=1)
  
#  l <-  applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
#                   as.formula("validate_true ~ (broad_mutect + dkfz + sanger)*(wgs_nvaf + wgs_tvaf + varlen + gencode + cosmic + dbsnp)"), 
#                   "stacked_logistic_regression", seed=1) 

  l$snvs$logisticRegression <- binarize(l$snvs$logistic_regression)
  l$snv_calls$logisticRegression <- binarize(l$snv_calls$logistic_regression)
  
  l$snvs$decisionTree <- binarize(l$snvs$decision_tree)
  l$snv_calls$decisionTree <- binarize(l$snv_calls$decision_tree)
  
#  l$snvs$stackedLogisticRegression <- binarize(l$snvs$stacked_logistic_regression, .4)
#  l$snv_calls$stackedLogisticRegression <- binarize(l$snv_calls$stacked_logistic_regression, .4)
  
  return(l)
}

simple_snv_models <- function(validated.calls, all.calls, seed=1) {
  l <-  applymodel(validated.calls, all.calls, glmlearn, glmpredict, 
                   as.formula("validate_true ~ wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_tvar_avgbaseposn + wgs_tvar_avgbaseq + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger + muse_feature"), 
                   "logistic_regression", seed=1)
  
  l <-  applymodel(l$snvs, l$snv_calls, ctree, ctreepredict, 
                   as.formula("validate_true ~ wgs_nvaf + wgs_tvaf + wgs_nvardepth + wgs_tvardepth + wgs_tvar_avgbaseposn + wgs_tvar_avgbaseq + gencode + cosmic + dbsnp + broad_mutect + dkfz + sanger + muse_feature"), 
                   "decision_tree", seed=1)
  
  l <-  applymodel(l$snvs, l$snv_calls, glmnetlearn, glmnetpredict, 
                   as.formula("validate_true ~ (broad_mutect + dkfz + sanger)*(wgs_nvaf + wgs_tvaf + gencode + cosmic + dbsnp + muse_feature)"), 
                   "stacked_logistic_regression", seed=1) 
  
  l$snvs$logisticRegression <- binarize(l$snvs$logistic_regression)
  l$snv_calls$logisticRegression <- binarize(l$snv_calls$logistic_regression)
  
  l$snvs$decisionTree <- binarize(l$snvs$decision_tree)
  l$snv_calls$decisionTree <- binarize(l$snv_calls$decision_tree)
  
  l$snvs$stackedLogisticRegression <- binarize(l$snvs$stacked_logistic_regression, .4)
  l$snv_calls$stackedLogisticRegression <- binarize(l$snv_calls$stacked_logistic_regression, .4)
  
  return(l)
}
doboxplot <- function(l) {
  boxplots(l$snvs, l$snv_calls, 
           callers=c('broad_mutect','dkfz', 'sanger', 'union', 'intersect2', 'intersect3', 'two_plus', 'logisticRegression', 'decisionTree')) 
}

dovafplot <- function(l, variant="SNV") {
  results <- corrected.accuracies.by.caller.by.vaf(l$snvs, l$snv_calls, 
           callers=c('union', 'intersect2', 'intersect3', 'two_plus', 'logisticRegression', 'decisionTree')) 
  results <- melt(results, id=c('caller','VAF'))
  results$feature <- FALSE
  results$feature[results$caller %in% c("logisticRegression","decisionTree","stackedLogisticRegression")] <- TRUE
  ggplot(results, aes(x=VAF,y=value)) + 
    geom_line(aes(group=caller,color=caller,linetype=feature)) + 
    facet_grid(variable~.)  +
    ggtitle(paste0("Accuracies of ", variant, " merge by VAF")) +
    ylab("Accuracy")
}

