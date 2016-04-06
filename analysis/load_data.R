#
# load in a CSV, and add useful columns such as:
#  - concordance
#  - 'seenby_' columns for each caller
#  - common_sample - if seen by all callers
#  - union/intersect2/intersect3 : derived callers from the core pipelines
#  - binned values for some continuous featuers: wgs_tvaf, homopolymer count, indel size
#
ingest_csv <- function(filename, callers, keep.lowdepth=FALSE) {
  data <- read.csv(filename)
  if (!keep.lowdepth) {
    data <- data[data$status != "LOWDEPTH",]
    data$status <- factor(data$status)
  }
  
  caller_columns <- sapply(callers, function(x) grep(x, names(data)))
  data$concordance <- rowSums(data[,caller_columns])
  data <- data[data$concordance>0, ]
  if ((sum(complete.cases(data$wgs_tvaf))>0) && (max(data$wgs_tvaf, na.rm=TRUE) > 0))
    data$binned_wgs_tvaf <- cut(data$wgs_tvaf, c(0,.1,.2,.3,.5,1), include.lowest=TRUE)
  data$validate_true <- data$status == "PASS"
  data$ref <- as.character(data$ref)
  data$alt <- as.character(data$alt)
  data$mnv <- ifelse( nchar(data$ref)>1, 1, 0 )
  data$varlen <- nchar(data$alt) - nchar(data$ref)
  
  col_names <- names(data)
  core_callers <- as.vector(na.omit(c(grep("broad",col_names), match("dkfz",col_names), match("embl_delly", col_names), match("sanger",col_names))))
  ncore <- rowSums(data[,core_callers])
  
  data$common_sample <- rep(0, nrow(data))
  for (caller in callers) {
    t <- tapply(data[[caller]], data$sample, FUN=sum)
    missed.samples <- dimnames(t)[[1]][t==0]
    seen <- paste0("seenby_",caller)
    data[[seen]] <- rep(1, nrow(data))
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
  
  data$repeat_masker[is.na(data$repeat_masker)] <- 0
  data$cosmic[is.na(data$cosmic)] <- 0
  data$dbsnp[is.na(data$dbsnp)] <- 0
  data$muse_feature[is.na(data$muse_feature)] <- 0
  
  data$indelsize <- abs(nchar(data$ref) - nchar(data$alt))
  if (max(data$indelsize, na.RM=TRUE) > 1)
    data$binned_indelsize <- cut(data$indelsize, c(0,3,5,10,25,50,100,250,Inf), include.lowest=TRUE, ordered_result=TRUE)
  if (sum(complete.cases(data$repeat_count))>0 && max(data$repeat_count, na.RM=TRUE) > 1)
    data$binned_homopolymer <- cut( data$repeat_count, c(0, 3, 10, 30, Inf), include.lowest=TRUE, ordered_result=TRUE)
  return(data)
}


array4_indel_calls <- ingest_csv('csvs/newmasters/array4_allcalls_indel.csv', indel_callers)
array3_indel_calls <- ingest_csv('csvs/newmasters/array3_allcalls_indel.csv', indel_callers)
array2_indel_calls <- ingest_csv('csvs/newmasters/array2_allcalls_indel.csv', indel_callers)
array1_indel_calls <- ingest_csv('csvs/newmasters/array1_allcalls_indel.csv', indel_callers)

array4_snv_calls <- ingest_csv('csvs/newmasters/array4_allcalls_snv_mnv.csv', snv_callers)
array3_snv_calls <- ingest_csv('csvs/newmasters/array3_allcalls_snv_mnv.csv', snv_callers)
array2_snv_calls <- ingest_csv('csvs/newmasters/array2_allcalls_snv_mnv.csv', snv_callers)
array1_snv_calls <- ingest_csv('csvs/newmasters/array1_allcalls_snv_mnv.csv', snv_callers)

#array4_sv_calls <- ingest_csv('csvs/newmasters/array4_allcalls_sv.csv', sv_callers)
#array3_sv_calls <- ingest_csv('csvs/newmasters/array3_allcalls_sv.csv', sv_callers)
#array2_sv_calls <- ingest_csv('csvs/newmasters/array2_allcalls_sv.csv', sv_callers)
#array1_sv_calls <- ingest_csv('csvs/newmasters/array1_allcalls_sv.csv', sv_callers)

array1_indels <- ingest_csv('csvs/newmasters/array1_indel.csv', indel_callers)
array2_indels <- ingest_csv('csvs/newmasters/array2_indel.csv', indel_callers)
array3_indels <- ingest_csv('csvs/newmasters/array3_indel.csv', indel_callers)
array4_indels <- ingest_csv('csvs/newmasters/array4_indel.csv', indel_callers)

array1_snvs <- ingest_csv('csvs/newmasters/array1_snv_mnv.csv', snv_callers)
array2_snvs <- ingest_csv('csvs/newmasters/array2_snv_mnv.csv', snv_callers)
array3_snvs <- ingest_csv('csvs/newmasters/array3_snv_mnv.csv', snv_callers)
array4_snvs <- ingest_csv('csvs/newmasters/array4_snv_mnv.csv', snv_callers)

#array4_svs <- ingest_csv('csvs/newmasters/array4_sv.csv', sv_callers)
#array3_svs <- ingest_csv('csvs/newmasters/array3_sv.csv', sv_callers)
#array2_svs <- ingest_csv('csvs/newmasters/array2_sv.csv', sv_callers)
#array1_svs <- ingest_csv('csvs/newmasters/array1_sv.csv', sv_callers)

snv_calls <- rbind(array1_snv_calls, array2_snv_calls, array3_snv_calls, array4_snv_calls)
indel_calls <- rbind(array1_indel_calls, array2_indel_calls, array3_indel_calls, array4_indel_calls)
#sv_calls <- rbind(array1_sv_calls, array2_sv_calls, array3_sv_calls, array4_sv_calls)

snvs <- rbind(array1_snvs, array2_snvs, array3_snvs, array4_snvs)
snvs$varlen <- 0
snvs$repeat_count <- 0

indels <- rbind(array1_indels, array2_indels, array3_indels, array4_indels)
indels$wgs_tvar_avgbaseq <- 0
indels$wgs_tvar_avgbaseposn <- 0
indel_calls$wgs_tvar_avgbaseq <- 0
indel_calls$wgs_tvar_avgbaseposn <- 0

#svs <- rbind(array1_svs, array2_svs, array3_svs, array4_svs)

snvs <- snvs[snvs$common_sample & snvs$repeat_masker == 0, ]
indels <- indels[indels$common_sample & indels$repeat_masker == 0, ]
snv_calls <- snv_calls[snv_calls$common_sample & snv_calls$repeat_masker == 0, ]
indel_calls <- indel_calls[indel_calls$common_sample & indel_calls$repeat_masker == 0, ]

snvs$gencode <- as.factor(ifelse(snvs$gencode=="", "", "gene"))
indels$gencode <- as.factor(ifelse(indels$gencode=="", "", "gene"))
snv_calls$gencode <- as.factor(ifelse(snv_calls$gencode=="", "", "gene"))
indel_calls$gencode <- as.factor(ifelse(indel_calls$gencode=="", "", "gene"))

snv_core_calls <- snv_calls[snv_calls$union == 1, ]
snvs_core <- snvs[snvs$union == 1, ]
indel_core_calls <- indel_calls[indel_calls$union == 1, ]
indels_core <- indels[indels$union == 1, ]

indel_norepeats <- indels[indels$repeat_count < 5,]
indel_calls_norepeats <- indel_calls[indel_calls$repeat_count < 5,]
indel_core_calls_norepeats <- indel_core_calls[indel_core_calls$repeat_count < 5,]

snv_core_calls$two_plus <- ifelse(snv_core_calls$broad_mutect + snv_core_calls$dkfz + snv_core_calls$sanger + snv_core_calls$muse_feature >= 2, 1, 0)
snv_calls$two_plus <- ifelse(snv_calls$broad_mutect + snv_calls$dkfz + snv_calls$sanger + snv_calls$muse_feature >= 2, 1, 0)
indel_core_calls$two_plus <- indel_core_calls$intersect2
indel_calls$two_plus <- indel_calls$intersect2
snvs_core$two_plus <- ifelse(snvs_core$broad_mutect + snvs_core$dkfz + snvs_core$sanger + snvs_core$muse_feature >= 2, 1, 0)
indels_core$two_plus <- indels_core$intersect2
snvs$two_plus <- ifelse(snvs$broad_mutect + snvs$dkfz + snvs$sanger + snvs$muse_feature >= 2, 1, 0)
indels$two_plus <- indels$intersect2

#rocplot(indels, indel_calls, formulae, names, indel_callers_plus_derived, indel_derived, "Indel Calls: Array 2+3+4: Corrected Accuracies, All")
#rocplot(indels_norepeats, indel_calls_norepeats, formulae, names, indel_callers_plus_derived, indel_derived, "Indel Calls: Array 2+3+4: Corrected Accuracies, repeat_count < 5")
#rocplot(snvs, snv_calls, formulae, names, snv_callers_plus_derived, snv_derived, "SNV Calls: Array 2+3+4: Corrected Accuracies, All")