snv_data <- function(filename) {
  snvs <- read.csv(filename)
  last_name_before_callers <- "indel_dist"
  
  firstcaller <- match(last_name_before_callers, names(snvs))+1
  lastcaller <- ncol(snvs)
  col_names <- names(snvs)
  callers <- col_names[firstcaller:lastcaller]
  
  for (caller in callers) {
    t <- tapply(snvs[[caller]], snvs$sample, FUN=sum)
    missed.samples <- dimnames(t)[[1]][t==0]
    seen <- paste0("seenby_",caller)
    snvs[[seen]] <- 1
    snvs[[seen]][snvs$sample %in% missed.samples] <- 0
  }
  col_names <- names(snvs)
  
  snvs$concordance <- rowSums(snvs[,firstcaller:lastcaller])
  if (max(snvs$wgs_tvaf, na.rm=TRUE) > 0)
    snvs$binned_wgs_tvaf <- cut(snvs$wgs_tvaf, c(0,.1,.2,.3,.5,1), include.lowest=TRUE)
  snvs$validate_true <- snvs$status == "PASS"
  snvs$ref <- as.character(snvs$ref)
  snvs$alt <- as.character(snvs$alt)
  snvs$mnv <- ifelse( nchar(snvs$ref)>1, 1, 0 )
  
  core_callers <- as.vector(na.omit(c(grep("broad",col_names), match("dkfz",col_names), match("embl_delly", col_names), match("sanger",col_names))))
  seenby_core_callers <- as.vector(na.omit(c(grep("seenby_broad",col_names), match("seenby_dkfz",col_names), match("seenby_embl_delly", col_names), match("seenby_sanger",col_names))))
  ncore <- rowSums(snvs[,core_callers])
  ncore_seenby <- rowSums(snvs[,seenby_core_callers])
  
  snvs$union <- ifelse( ncore > 0, 1, 0)
  snvs$intersect2 <- ifelse( ncore > 1, 1, 0)
  snvs$intersect3 <- ifelse( ncore > 2, 1, 0)
  
  snvs$indelsize <- abs(nchar(snvs$ref) - nchar(snvs$alt))
  if (max(snvs$indelsize, na.RM=TRUE) > 1)
    snvs$binned_indelsize <- cut(snvs$indelsize, c(0,3,5,10,25,50,100,250,max(snvs$indelsize)), include.lowest=TRUE, ordered_result=TRUE)
  if (max(snvs$repeat_count, na.RM=TRUE) > 1)
    snvs$binned_homopolymer <- cut( snvs$repeat_count, c(0, 3, 10, max(snvs$repeat_count)), include.lowest=TRUE, ordered_result=TRUE)
  return(snvs)
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

  if(is.factor(seendata$variable)) {
    vals <- levels(seendata$variable)
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