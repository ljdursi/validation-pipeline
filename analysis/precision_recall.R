derived <- c("union", "intersect2", "intersect3")

snv_callers <- c("adiscan", "broad_mutect", "dkfz", "lohcomplete", "mda_hgsc_gatk_muse", "oicr_bl", "oicr_sga", "sanger", "smufin", "wustl")
snv_callers_plus_derived <- c(snv_callers, derived)
snv_derived <- c(rep(FALSE, length(snv_callers)), rep(TRUE, length(derived)))

indel_callers <- c("broad_mutect", "crg_clindel", "dkfz", "novobreak", "oicr_sga", "sanger", "smufin", "wustl")
indel_callers_plus_derived <- c(indel_callers, derived)
indel_derived <- c(rep(FALSE, length(indel_callers)), rep(TRUE, length(derived)))

sv_callers <- c("broad_merged", "destruct", "embl_delly", "novobreak", "oicr_bl", "sanger", "smufin", "wustl")
sv_callers_plus_derived <- c(sv_callers, derived)
sv_derived <- c(rep(FALSE, length(sv_callers)), rep(TRUE, length(derived)))

core_callers_formula <- "validate_true ~ broad_mutect + dkfz + sanger + wgs_tvaf + wgs_nvaf"


#
# raw accuracies, not taking into account the selection procedure
#
accuracies <- function(validate_true, calls, weights=rep(1,length(calls))) {
  #ntruepos <- sum((validate_true==TRUE & calls==1)*weights)
  #nfalsepos <- sum((validate_true==FALSE & calls==1)*weights)
  #nfalseneg <- sum((validate_true==TRUE & calls==0)*weights)
  
  ntruepos <- sum(weights[validate_true==TRUE & calls==1])
  nfalsepos <- sum(weights[validate_true==FALSE & calls==1])
  nfalseneg <- sum(weights[validate_true==TRUE & calls==0])
  
  print(paste("ntruepos = ", ntruepos))
  print(paste("nfalsepos = ", nfalsepos))
  print(paste("nfalseneg = ", nfalseneg))
  sensitivity <- ntruepos/(ntruepos+nfalseneg)
  precision <- ntruepos/(ntruepos+nfalsepos)
  f1 <- 2*ntruepos/(2*ntruepos+nfalsepos+nfalseneg)
  return(list(sensitivity=sensitivity, precision=precision, f1=f1))
}

raw_accuracies_by_var <- function(data, caller, variable=NULL) {
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

uncorrected.accuracies <- function(validated.calls, all.calls, caller) {
  ntruepos <- sum(validated.calls$validate_true==TRUE & validated.calls[[caller]]==1)
  nfalsepos <- sum(validated.calls$validate_true==FALSE & validated.calls[[caller]]==1)
  nfalseneg <- sum(validated.calls$validate_true==TRUE & validated.calls[[caller]]==0)
  
  sensitivity <- ntruepos/(ntruepos+nfalseneg)
  precision <- ntruepos/(ntruepos+nfalsepos)
  f1 = 2*sensitivity*precision/(sensitivity+precision)
  data.frame(sensitivity=sensitivity, precision=precision, f1=f1, caller=caller)
}

corrected.accuracies <- function(validated.calls, all.calls, caller, 
                                 by=c("concordance","sample"), 
                                 combine=c("concordance","sample")) {

  validated.calls$concordance <- factor(validated.calls$concordance)
  validated.calls$sample <- factor(validated.calls$sample)
  
  all.calls <- all.calls[all.calls$sample %in% validated.calls$sample, ]
  all.calls$sample <- factor(all.calls$sample)
  all.calls$concordance <- factor(all.calls$concordance)
  
  form <- paste("~",paste(by, collapse='+'))
  
  all.positive.calls <- as.matrix(xtabs(as.formula(form), data=all.calls, subset=all.calls[[caller]]==1, drop.unused.levels=FALSE))
  all.negative.calls <- as.matrix(xtabs(as.formula(form), data=all.calls, subset=all.calls[[caller]]==0, drop.unused.levels=FALSE))
  
  validated.negative.calls <- as.matrix(xtabs(as.formula(form), data=validated.calls, subset=validated.calls[[caller]]==0, drop.unused.levels=FALSE))
  validated.tp <- as.matrix(xtabs(as.formula(form), data=validated.calls, subset=(validated.calls[[caller]]==1 & validated.calls$validate_true), drop.unused.levels=FALSE))
  validated.fn <- as.matrix(xtabs(as.formula(form), data=validated.calls, subset=(validated.calls[[caller]]==0 & validated.calls$validate_true), drop.unused.levels=FALSE))
  validated.fp <- as.matrix(xtabs(as.formula(form), data=validated.calls, subset=(validated.calls[[caller]]==1 & !validated.calls$validate_true), drop.unused.levels=FALSE))
  
  validated.positive.calls <- validated.tp + validated.fp

  stopifnot(dim(all.positive.calls)==dim(all.negative.calls))
  stopifnot(dim(all.positive.calls)==dim(validated.positive.calls))
  stopifnot(dim(all.positive.calls)==dim(validated.negative.calls))
  
  combined.dims <- sapply(combine, function(x) grep(x,by))
  
  eps <- 1.e-5
  pos.weights <- all.positive.calls/(validated.positive.calls + eps)
  neg.weights <- all.negative.calls/(validated.negative.calls + eps)
  
  if (length(combined.dims) == length(dim(pos.weights))) {
    marginsum <- function(x) sum(x)
  } else if (combined.dims[1] == 1) {
    marginsum <- function(x) colSums(x)
  } else {
    marginsum <- function(x) rowSums(x)
  }
  
  estimated.all.tp <- marginsum(validated.tp * pos.weights)
  estimated.all.fp <- marginsum(validated.fp * pos.weights)
  estimated.all.fn <- marginsum(validated.fn * neg.weights)

  sensitivity <- estimated.all.tp/(estimated.all.tp + estimated.all.fn+eps)
  precision <- estimated.all.tp/(estimated.all.tp + estimated.all.fp+eps)
  f1 <- 2*sensitivity*precision/(sensitivity+precision+eps)
  df <- data.frame(sensitivity=sensitivity, precision=precision, f1=f1, caller=caller)
  if (length(combine) != length(by)) {
    if (!("sample" %in% combine)) {
      df$sample <- names(estimated.all.tp)
    } else {
      df$concordance <- names(estimated.all.tp)
    }
  }
  return(df)
}

corrected.accuracies.by.caller <- function(validated.call.data, all.call.data, callers) {
  results <- do.call(rbind, lapply(callers, function(x) corrected.accuracies(validated.call.data, all.call.data, x)))
  return(results)
}

uncorrected.accuracies.by.caller <- function(validated.call.data, all.call.data, callers) {
  results <- do.call(rbind, lapply(callers, function(x) uncorrected.accuracies(validated.call.data, all.call.data, x)))
  return(results)
}

corrected.accuracies.by.caller.pooled <- function(validated.call.data, all.call.data, callers) {
  results <- do.call(rbind, lapply(callers, function(x) corrected.accuracies(validated.call.data, all.call.data, x, by=c("concordance"), combine=c("concordance"))))
  return(results)
}

corrected.accuracies.by.caller.by.sample <- function(validated.call.data, all.call.data, callers) {
  results <- do.call(rbind, lapply(callers, function(x) corrected.accuracies(validated.call.data, all.call.data, x, by=c("concordance","sample"), combine=c("concordance"))))
  return(results)
}

corrected.accuracies.by.caller.by.vaf <- function(validated.call.data, all.call.data, callers) {
  results <- do.call(rbind, lapply(callers, function(x) corrected.accuracies(validated.call.data, all.call.data, x, by=c("concordance","binned_wgs_tvaf"), combine=c("concordance"))))
  colnames(results) <- c("sensitivity","precision","f1","caller","VAF")
  results$VAF <- factor(results$VAF, levels=levels(validated.call.data$binned_wgs_tvaf))
  return(results)
}

corrected.accuracies.by.caller.by.indelsize <- function(validated.call.data, all.call.data, callers) {
  results <- do.call(rbind, lapply(callers, function(x) corrected.accuracies(validated.call.data, all.call.data, x, by=c("concordance","binned_indelsize"), combine=c("concordance"))))
  colnames(results) <- c("sensitivity","precision","f1","caller","indelsize")
  results$indelsize <- factor(results$indelsize, levels=levels(validated.call.data$binned_indelsize))
  return(results)
}


corrected.accuracies.with.cis.by.caller <- function(validated.call.data, all.call.data, callers, ntrials=100) {
  validated.calls <- validated.call.data[,c("concordance","validate_true", "sample", callers)]
  n <- nrow(validated.calls)
  
  results <- do.call(rbind, lapply(1:ntrials, function(x) {
                                                samp.validated <- validated.calls[sample.int(n, replace=TRUE),];
                                                df <- corrected.accuracies.by.caller(samp.validated, all.call.data, callers);
                                                df$trial <- x;
                                                return(df);
                                              }))
  means <- aggregate(results[,c("sensitivity","precision","f1")], by=list(results$caller), FUN=mean)
  colnames(means)[1] <- "caller"
   
  lowci <- aggregate(results[,c("sensitivity","precision","f1")], by=list(results$caller), FUN=function(x) quantile(x,probs=c(0.05)))                  
  colnames(lowci) <- c("caller", "sensitivity_low_ci", "precision_low_ci", "f1_low_ci")
  
  hici <- aggregate(results[,c("sensitivity","precision","f1")], by=list(results$caller), FUN=function(x) quantile(x,probs=c(0.95)))                  
  colnames(hici) <- c("caller", "sensitivity_hi_ci", "precision_hi_ci", "f1_hi_ci")
  
  results <- merge(means, merge(lowci, hici))
  return(results)
}

printcis <- function(data.with.cis) {
  data.with.cis$caller <- as.character(data.with.cis$caller)
  row_to_string <- function(x) {
    sprintf("%20s: sensitivity = %5.3f (%8.5f-%8.5f)\n                    : precision   = %5.3f (%8.5f-%8.5f)\n                    : f1          = %5.3f (%8.5f-%8.5f)\n", 
            x["caller"], as.numeric(x["sensitivity"]), as.numeric(x["sensitivity_low_ci"]), as.numeric(x["sensitivity_hi_ci"]), 
            as.numeric(x["precision"]), as.numeric(x["precision_low_ci"]), as.numeric(x["precision_hi_ci"]), 
            as.numeric(x["f1"]), as.numeric(x["f1_low_ci"]), as.numeric(x["f1_hi_ci"]))
  }
  as.vector(apply(data.with.cis, 1, row_to_string))
}
