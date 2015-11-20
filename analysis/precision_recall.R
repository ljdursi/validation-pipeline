snv_callers <- c("adiscan", "broad_mutect", "dkfz", "lohcomplete", "mda_hgsc_gatk_muse", "oicr_bl", "oicr_sga", "sanger", "smufin", "wustl")
indel_callers <- c("broad_mutect", "crg_clindel", "dkfz", "novobreak", "oicr_sga", "sanger", "smufin", "wustl")
sv_callers <- c("broad_merged", "destruct", "embl_delly", "novobreak", "sanger", "smufin")

csv_data <- function(filename, callers) {
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
  
  for (caller in callers) {
    t <- tapply(data[[caller]], data$sample, FUN=sum)
    missed.samples <- dimnames(t)[[1]][t==0]
    seen <- paste0("seenby_",caller)
    data[[seen]] <- 1
    data[[seen]][data$sample %in% missed.samples] <- 0
  }
  col_names <- names(data)
  seenby_core_callers <- as.vector(na.omit(c(grep("seenby_broad",col_names), match("seenby_dkfz",col_names), match("seenby_embl_delly", col_names), match("seenby_sanger",col_names))))
  ncore_seenby <- rowSums(data[,seenby_core_callers])
  
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

accuracies <- function(validate_true, calls) {
  ntruepos = sum(validate_true==TRUE & calls==1)
  nfalsepos = sum(validate_true==FALSE & calls==1)
  nfalseneg = sum(validate_true==TRUE & calls==0)
  sensitivity = ntruepos/(ntruepos+nfalseneg)
  precision = ntruepos/(ntruepos+nfalsepos)
  f1 = 2*ntruepos/(2*ntruepos+nfalsepos+nfalseneg)
  return(list(sensitivity=sensitivity, precision=precision, f1=f1))
}

accuracies_by_var <- function(data, caller, variable=NULL) {
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

caller_sensitivities <- function(data, callers, variable="binned_wgs_tvaf") {
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
rocplot <- function(data, formulae, names, callers, derived, title) {
  s <- caller_sensitivities(data, callers, NULL)
  s$derived <- derived
  
  allmodels <- data.frame()
  for (i in 1:length(names)) {
    for (trial in 1:10) {
    l <- split(data, sample(1:2, nrow(data), replace=TRUE))
    test <- l[[1]]; train <- l[[2]]
    
    treemodel <- ctree(as.formula(formulae[i]), data=train)
    results <- modelROC(test, treemodel)
    results$model <- names[i]
    results$type <- "tree"
    
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

table.to.vect <- function(values, maxcols=12) {
  vect <- rep(0, maxcols)
  t <- table(values)
  vect[as.integer(names(t))] <- t
  return(vect)
}

corrected.accuracies.per.sample <- function(selected.call.data, all.call.data, callername, samplename) {
  verbose <- FALSE
  selected.calls <- selected.call.data[selected.call.data$sample==samplename, ]
  all.calls <- all.call.data[all.call.data$sample==samplename,]
  
  all.positive.calls.by.concordance <- table.to.vect(all.calls[all.calls[[callername]]==1,]$concordance)
  all.negative.calls.by.concordance <- table.to.vect(all.calls[all.calls[[callername]]==0,]$concordance)
  
  if (verbose) {
    print("all.positive.calls.by.concordance")
    print(all.positive.calls.by.concordance)
    print("all.negative.calls.by.concordance")
    print(all.negative.calls.by.concordance)
  }
  
  validated.tp.by.concordance <- table.to.vect(selected.calls[selected.calls[[callername]]==1 & selected.calls$validate_true,]$concordance)
  validated.fp.by.concordance <- table.to.vect(selected.calls[selected.calls[[callername]]==1 & !selected.calls$validate_true,]$concordance)
  validated.tn.by.concordance <- table.to.vect(selected.calls[selected.calls[[callername]]==0 & !selected.calls$validate_true,]$concordance)
  validated.fn.by.concordance <- table.to.vect(selected.calls[selected.calls[[callername]]==0 & selected.calls$validate_true,]$concordance)

  if (verbose) {
    print("validated.tp.by.concordance")
    print(validated.tp.by.concordance)
    print("validated.fp.by.concordance")
    print(validated.fp.by.concordance)
    print("validated.tn.by.concordance")
    print(validated.tn.by.concordance)
    print("validated.fn.by.concordance")
    print(validated.fn.by.concordance)
  }
  validated.positive.calls.by.concordance <- validated.tp.by.concordance + validated.fp.by.concordance
  validated.negative.calls.by.concordance <- validated.tn.by.concordance + validated.fn.by.concordance
  
  if (verbose) {
    print("validated positive calls")
    print(validated.positive.calls.by.concordance)
    print("validated negative calls")
    print(validated.negative.calls.by.concordance)
  }
  eps <- 0.001 #(to avoid getting NaNs when we have 0/0; we want those terms to not contribute, so give 0)
  
  if (verbose) print("estimated all tp")
  estimated.all.tp <- validated.tp.by.concordance/(validated.positive.calls.by.concordance+eps) * all.positive.calls.by.concordance
  if (verbose) print(estimated.all.tp)
  estimated.all.tp <- sum(estimated.all.tp)
  if (verbose) print(estimated.all.tp)
  
  if (verbose) print("estimated all fp")
  estimated.all.fp <- validated.fp.by.concordance/(validated.positive.calls.by.concordance+eps) * all.positive.calls.by.concordance
  if (verbose) print(estimated.all.fp)
  estimated.all.fp <- sum(estimated.all.fp)
  if (verbose) print(estimated.all.fp)
  
  if (verbose) print("estimated all fn")
  estimated.all.fn <- validated.fn.by.concordance/(validated.negative.calls.by.concordance+eps) * all.negative.calls.by.concordance
  if (verbose) print(estimated.all.fn)
  estimated.all.fn <- sum(estimated.all.fn)
  if (verbose) print(estimated.all.fn)
  
  sensitivity <- estimated.all.tp/(estimated.all.tp + estimated.all.fn)
  precision <- estimated.all.tp/(estimated.all.tp + estimated.all.fp)
  
  f1 = 2*sensitivity*precision/(sensitivity+precision)
  return(list(sensitivity=sensitivity, precision=precision, f1=f1))
}

corrected.accuracies.by.sample <- function(selected.call.data, all.call.data, callername) {
  vals <- levels(selected.call.data$sample)
  results <- lapply(vals, 
                    function(x) corrected.accuracies.per.sample(selected.call.data, all.call.data, callername, x))
  df <- data.frame(t(matrix(unlist(results), ncol=length(vals))))
  colnames(df) <- c("sensitivity","precision","f1")
  df$sample <- vals
  df$caller <- callername
  return(df)
}