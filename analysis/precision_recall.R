snv_data <- function(filename) {
  snvs <- read.csv(filename)
  snvs$concordance <- snvs$adiscan + snvs$broad_mutect + snvs$dkfz + snvs$lohcomplete + snvs$mda_hgsc_gatk_muse + snvs$oicr_bl + snvs$oicr_sga + snvs$sanger + snvs$smufin + snvs$wustl
  snvs$binned_wgs_tvaf <- cut(snvs$wgs_tvaf, c(0,.1,.2,.3,.5,1), include.lowest=TRUE)
  snvs$validate_true <- snvs$status == "PASS"
  snvs$ref <- as.character(snvs$ref)
  snvs$alt <- as.character(snvs$alt)
  snvs$mnv <- ifelse( nchar(snvs$ref)>1, 1, 0 )
  snvs$union <- ifelse( snvs$broad_mutect + snvs$dkfz + snvs$sanger > 0, 1, 0)
  snvs$intersect2 <- ifelse( snvs$broad_mutect + snvs$dkfz + snvs$sanger > 1, 1, 0)
  snvs$intersect3 <- ifelse( snvs$broad_mutect + snvs$dkfz + snvs$sanger > 2, 1, 0)
  return(snvs)
}

accuracies <- function(data, caller) {
  ntruepos = sum(data$validate_true==1 & data[[caller]]==1)
  nfalsepos = sum(data$validate_true==0 & data[[caller]]==1)
  nfalseneg = sum(data$validate_true==1 & data[[caller]]==0)
  sensitivity = ntruepos/(ntruepos+nfalseneg)
  precision = ntruepos/(ntruepos+nfalsepos)
  f1 = 2*ntruepos/(2*ntruepos+nfalsepos+nfalseneg)
  return(list(sensitivity=sensitivity, precision=precision, f1=f1))
}

accuracies_by_var <- function(data, caller, variable=NULL) {
  if (is.null(variable)) 
    return( list(all=accuracies(data,caller)) )
  
  vals <- unique(data[[variable]])
  valnames <- as.character((vals))
  
  results <- sapply(unique(data[[variable]]), function(val) accuracies(data[data[[variable]]==val,],caller))
  colnames(results) <- valnames
  return(results)
}

caller_sensitivities <- function(filename) {
  snvs <- snv_data(filename)
  callers <- names(snvs)[11:20]
  for (caller in callers) {
    print(caller)
    r <- accuracies_by_var(snvs, caller, "binned_wgs_tvaf")
    print(r)
  }
}

