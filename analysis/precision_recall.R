snv_data <- function(filename) {
  snvs <- read.csv(filename)
  last_name_before_callers <- "indel_dist"
  
  firstcaller <- match(last_name_before_callers, names(snvs))+1
  lastcaller <- ncol(snvs)
  col_names <- names(snvs)
  callers <- col_names[firstcaller:lastcaller]
  
  snvs$concordance <- rowSums(snvs[,firstcaller:lastcaller])
  #snvs$binned_wgs_tvaf <- cut(snvs$wgs_tvaf, c(0,.1,.2,.3,.5,1), include.lowest=TRUE)
  snvs$validate_true <- snvs$status == "PASS"
  snvs$ref <- as.character(snvs$ref)
  snvs$alt <- as.character(snvs$alt)
  snvs$mnv <- ifelse( nchar(snvs$ref)>1, 1, 0 )
  
  core_callers <- na.omit(c(grep("broad",col_names), match("dkfz",col_names), match("embl_delly", col_names), match("sanger",col_names)))
  ncore <- rowSums(snvs[,core_callers])
  
  snvs$union <- ifelse( ncore > 0, 1, 0)
  snvs$intersect2 <- ifelse( ncore > 1, 1, 0)
  snvs$intersect3 <- ifelse( ncore  > 2, 1, 0)
  
  snvs$indelsize <- abs(nchar(snvs$ref) - nchar(snvs$alt))
  #snvs$binned_indelsize <- cut(snvs$indelsize, c(0,3,5,10,25,50,100,250,max(snvs$indelsize)), include.lowest=TRUE, ordered_result=TRUE)
  #snvs$binned_homopolymer <- cut( snvs$repeat_count, c(0, 3, 10, max(snvs$repeat_count)), include.lowest=TRUE, ordered_result=TRUE)
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
  if (is.null(variable)) {
    df <- data.frame( accuracies(data$validate_true, data[[caller]]) )
    df$caller <- caller
    return(df)
  }

  if(is.factor(data$variable)) {
    vals <- levels(data$variable)
  } else {
    vals <- sort(unique(data[[variable]]))
  }
  valnames <- as.character((vals))
  
  results <- lapply(vals, function(v) { rows <- data[[variable]]==v;
                                        accuracies(data[rows,]$validate_true,
                                                   data[[caller]][rows]) } ) 
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

library(party)
library(randomForest)
rocplot <- function(data, formulae, names, callers, derived) {
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
    }
  }
  ggplot() + geom_point(data=s, aes(x=sensitivity, y=precision, label=caller)) + 
             geom_text(data=s, aes(x=sensitivity, y=precision, label=caller, angle=-45), hjust=0, vjust=0) + 
             geom_smooth(data=allmodels, aes(x=sensitivity, y=precision, color=type, linetype=model))
}