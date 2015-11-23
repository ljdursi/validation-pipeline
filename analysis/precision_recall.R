derived <- c("union", "intersect2", "intersect3")

snv_callers <- c("adiscan", "broad_mutect", "dkfz", "lohcomplete", "mda_hgsc_gatk_muse", "oicr_bl", "oicr_sga", "sanger", "smufin", "wustl")
snv_callers_plus_derived <- c(snv_callers, derived)
snv_derived <- c(rep(FALSE, length(snv_callers)), rep(TRUE, length(derived)))
indel_callers <- c("broad_mutect", "crg_clindel", "dkfz", "novobreak", "oicr_sga", "sanger", "smufin", "wustl")
indel_callers_plus_derived <- c(indel_callers, derived)
indel_derived <- c(rep(FALSE, length(indel_callers)), rep(TRUE, length(derived)))
sv_callers <- c("broad_merged", "destruct", "embl_delly", "novobreak", "sanger", "smufin")
sv_callers_plus_derived <- c(snv_callers, derived)
sv_derived <- c(rep(FALSE, length(sv_callers)), rep(TRUE, length(derived)))

core_callers_formula <- "validate_true ~ broad_mutect + dkfz + sanger + wgs_tvaf + wgs_nvaf"

#
# load in a CSV, and add useful columns such as:
#  - concordance
#  - 'seenby_' columns for each caller
#  - common_sample - if seen by all callers
#  - union/intersect2/intersect3 : derived callers from the core pipelines
#  - binned values for some continuous featuers: wgs_tvaf, homopolymer count, indel size
#
load_csv <- function(filename, callers) {
  data <- read.csv(filename)

  caller_columns <- sapply(callers, function(x) grep(x, names(data)))
  data$concordance <- rowSums(data[,caller_columns])
  if (max(data$wgs_tvaf, na.rm=TRUE) > 0)
    data$binned_wgs_tvaf <- cut(data$wgs_tvaf, c(0,.1,.2,.3,.5,1), include.lowest=TRUE)
  data$validate_true <- data$status == "PASS"
  data$ref <- as.character(data$ref)
  data$alt <- as.character(data$alt)
  data$mnv <- ifelse( nchar(data$ref)>1, 1, 0 )
  
  col_names <- names(data)
  core_callers <- as.vector(na.omit(c(grep("broad",col_names), match("dkfz",col_names), match("embl_delly", col_names), match("sanger",col_names))))
  ncore <- rowSums(data[,core_callers])
  
  data$common_sample <- 0
  for (caller in callers) {
    t <- tapply(data[[caller]], data$sample, FUN=sum)
    missed.samples <- dimnames(t)[[1]][t==0]
    seen <- paste0("seenby_",caller)
    data[[seen]] <- 1
    data[[seen]][data$sample %in% missed.samples] <- 0
    data$common_sample <- data$common_sample + data[[seen]]
  }
  col_names <- names(data)
  seenby_core_callers <- as.vector(na.omit(c(grep("seenby_broad",col_names), match("seenby_dkfz",col_names), match("seenby_embl_delly", col_names), match("seenby_sanger",col_names))))
  ncore_seenby <- rowSums(data[,seenby_core_callers])
  
  data$common_sample <- data$common_sample == length(callers)
  
  data$union <- ifelse( ncore > 0, 1, 0)
  data$intersect2 <- ifelse( ncore > 1, 1, 0)
  data$intersect3 <- ifelse( ncore > 2, 1, 0)
  
  data$indelsize <- abs(nchar(data$ref) - nchar(data$alt))
  if (max(data$indelsize, na.RM=TRUE) > 1)
    data$binned_indelsize <- cut(data$indelsize, c(0,3,5,10,25,50,100,250,max(data$indelsize)), include.lowest=TRUE, ordered_result=TRUE)
  if (max(data$repeat_count, na.RM=TRUE) > 1)
    data$binned_homopolymer <- cut( data$repeat_count, c(0, 3, 10, max(data$repeat_count)), include.lowest=TRUE, ordered_result=TRUE)
  return(data)
}

#
# raw accuracies, not taking into account the selection procedure
#
rawaccuracies <- function(validate_true, calls) {
  ntruepos = sum(validate_true==TRUE & calls==1)
  nfalsepos = sum(validate_true==FALSE & calls==1)
  nfalseneg = sum(validate_true==TRUE & calls==0)
  sensitivity = ntruepos/(ntruepos+nfalseneg)
  precision = ntruepos/(ntruepos+nfalsepos)
  f1 = 2*ntruepos/(2*ntruepos+nfalsepos+nfalseneg)
  return(list(sensitivity=sensitivity, precision=precision, f1=f1))
}

rawaccuracies_by_var <- function(data, caller, variable=NULL) {
  seenby_caller <- paste0("seenby_",caller)
  if (seenby_caller %in% names(data))
    seendata <- data[data[[seenby_caller]]==1,]
  else
    seendata <- data
  
  if (is.null(variable)) {
    df <- data.frame( accuracies(seendata$validate_true, seendata[[caller]]) )
    df$caller <- caller
    return(df)
  }

  if(is.factor(seendata[[variable]])) {
    vals <- levels(seendata[[variable]])
  } else {
    vals <- sort(unique(seendata[[variable]]))
  }
  valnames <- as.character((vals))
  

  results <- lapply(vals, function(v) { rows <- seendata[[variable]]==v;
                                        accuracies(seendata[rows,]$validate_true,
                                                   seendata[[caller]][rows]) } )

  df <- data.frame(t(matrix(unlist(results), ncol=length(vals))))
  colnames(df) <- c("sensitivity","precision","f1")
  df[[variable]] <- vals
  df$caller <- caller

  return(df)
}

raw_caller_sensitivities <- function(data, callers, variable="binned_wgs_tvaf") {
  do.call(rbind, lapply(callers, function(x) accuracies_by_var(data, x, variable)))
}

modelROC <- function(data, model) {
  vals <- (10:95)/100
  prediction <- predict(model, newdata=data)
  results <- sapply(vals, function(x) accuracies(data$validate_true, prediction > x))
  df <- data.frame(t(matrix(unlist(results),ncol=length(vals))))
  colnames(df) <- c("sensitivity","precision","f1")
  df$thresh <- vals
  return(df)
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
  s <- corrected.accuracies(data, allcalls, callers)
  s$derived <- derived
  
  allmodels <- data.frame()
  for (i in 1:length(names)) {
    for (trial in 1:10) {
    l <- split(data, sample(1:2, nrow(data), replace=TRUE))
    test <- l[[1]]; train <- l[[2]]
    
    treemodel <- ctree(as.formula(formulae[i]), data=train)
    results <- modelROC(test, treemodel)
    results$model <- names[i]
    results$type <- "decision tree"
    
    allmodels <- rbind(allmodels, results)
    
    glmmodel <- glm(as.formula(formulae[i]), data=train)
    results <- modelROC(test, glmmodel)
    results$model <- names[i]
    results$type <- "logistic regression"
    
    allmodels <- rbind(allmodels, results)
    
    rfmodel <- rf(formulae[i], train)
    results <- rfROC(test, rfmodel)
    results$model <- names[i]
    results$type <- "random forest"
    
    allmodels <- rbind(allmodels, results)
    }
  }
  
  f1s <- seq(.6,.95,.05)
  roc.contours <- const.f1.rocs(f1s)
  ggplot() + geom_point(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller)) +
             geom_point(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller), color='blue') +
             geom_text(data=s[s$derived==FALSE,], aes(x=sensitivity, y=precision, label=caller), angle=45, hjust=1, vjust=1) + 
             geom_text(data=s[s$derived==TRUE,], aes(x=sensitivity, y=precision, label=caller, angle=+45), hjust=1, vjust=1, color='blue') + 
             geom_smooth(data=allmodels, aes(x=sensitivity, y=precision, color=type, linetype=model)) +
             geom_line(data=roc.contours, aes(x=sensitivity, y=precision, group=f1), color='grey', alpha=0.3) +
             ggtitle(title)
  
}

allseen <- function(data) {
  allseen_cols <- grep("seenby_", names(data))
  numseens <- rowSums(data[,allseen_cols])
  return(data[numseens==max(numseens),])
}

estimated.num.mutations <- function(selected.call.data, all.call.data, samplename) {
  selected.calls <- selected.call.data[selected.call.data$sample==samplename,]
  all.calls <- all.call.data[all.call.data$sample==samplename,]
  
  total.calls.by.concordance <- table(all.calls$concordance)
  validated.calls.by.concordance <- table(selected.calls$concordance)
  validated.true.calls.by.concordance <- table(selected.calls$concordance[selected.calls$validate_true])
  
  truefrac.by.concordance <- validated.true.calls.by.concordance / validated.calls.by.concordance
  return(as.integer(sum(total.calls.by.concordance*truefrac.by.concordance)))
}

corrected.accuracies.by.caller <- function(validated.call.data, all.call.data, caller) {
  verbose <- FALSE
  validated.calls <- subset(validated.call.data, common_sample)
  validated.calls$concordance <- factor(validated.calls$concordance)
  validated.calls$sample <- factor(validated.calls$sample)
  commonsamples <- unique(validated.calls$sample)
  
  all.calls <- all.call.data[all.call.data$sample %in% commonsamples,]
  all.calls <- all.calls[all.calls$concordance > 0, ]
  all.calls$sample <- factor(all.calls$sample)
  all.calls$concordance <- factor(all.calls$concordance)
  
  all.positive.calls <- as.matrix(xtabs(~concordance+sample, data=all.calls, subset=all.calls[[caller]]==1, drop.unused.levels=FALSE))
  all.negative.calls <- as.matrix(xtabs(~concordance+sample, data=all.calls, subset=all.calls[[caller]]==0, drop.unused.levels=FALSE))
  
  validated.tp <- as.matrix(xtabs(~concordance+sample, data=validated.calls, subset=(validated.calls[[caller]]==1 & validated.calls$validate_true), drop.unused.levels=FALSE))
  validated.fp <- as.matrix(xtabs(~concordance+sample, data=validated.calls, subset=(validated.calls[[caller]]==1 & !validated.calls$validate_true), drop.unused.levels=FALSE))
  validated.tn <- as.matrix(xtabs(~concordance+sample, data=validated.calls, subset=(validated.calls[[caller]]==0 & !validated.calls$validate_true), drop.unused.levels=FALSE))
  validated.fn <- as.matrix(xtabs(~concordance+sample, data=validated.calls, subset=(validated.calls[[caller]]==0 & validated.calls$validate_true), drop.unused.levels=FALSE))

  stopifnot(dim(all.positive.calls)==dim(all.negative.calls))
  stopifnot(dim(all.positive.calls)==dim(validated.tp))
  stopifnot(dim(all.positive.calls)==dim(validated.fp))
  stopifnot(dim(all.positive.calls)==dim(validated.tn))
  stopifnot(dim(all.positive.calls)==dim(validated.fn))
  
  eps <- 0.0001
  estimated.all.tp <- validated.tp/(validated.tp+validated.fp+eps) * all.positive.calls
  estimated.all.fp <- validated.fp/(validated.tp+validated.fp+eps) * all.positive.calls
  estimated.all.fn <- validated.fn/(validated.tn+validated.fn+eps) * all.negative.calls
  estimated.all.tn <- validated.tn/(validated.tn+validated.fn+eps) * all.negative.calls
  
  #sensitivity.by.sample <- colSums(estimated.all.tp)/(colSums(estimated.all.tp + estimated.all.fn)+eps)
  #precision.by.sample <- colSums(estimated.all.tp)/(colSums(estimated.all.tp + estimated.all.fp)+eps)
  #f1.by.sample = 2*sensitivity.by.sample*precision.by.sample/(sensitivity.by.sample+precision.by.sample+eps)
  #print(sensitivity.by.sample)
  #print(precision.by.sample)
  #print(f1.by.sample)

  sensitivity <- sum(estimated.all.tp)/(sum(estimated.all.tp + estimated.all.fn)+eps)
  precision <- sum(estimated.all.tp)/(sum(estimated.all.tp + estimated.all.fp)+eps)
  f1 <- 2*sensitivity*precision/(sensitivity+precision+eps)
  return(list(sensitivity=sensitivity, precision=precision, f1=f1))
}

corrected.accuracies <- function(validated.call.data, all.call.data, callers) {
  results <- sapply(callers, function(x) corrected.accuracies.by.caller(validated.call.data, all.call.data, x))
  df <- data.frame(t(matrix(unlist(results), ncol=length(callers))))
  colnames(df) <- dimnames(results)[[1]]
  df$caller <- dimnames(results)[[2]]
  return(df)
}