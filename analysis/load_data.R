library(dplyr)

#
# load in a CSV, and add useful columns such as:
#  - concordance
#  - 'seenby_' columns for each caller
#  - common_sample - if seen by all callers
#  - union/intersect2/intersect3 : derived callers from the core pipelines
#  - binned values for some continuous featuers: wgs_tvaf, homopolymer count, indel size
#

find_columns <- function(names, desired_columns) {
  candidate_columns <- sapply(desired_columns, function(x) grep(x, names))
  cols <- c()
  for (col in candidate_columns) {
    cols<- c(cols, col)
  }
  return(cols)
}

ingest_csv <- function(filename, callers, keep.lowdepth=FALSE, core_callers=c('broad_mutect', 'dkfz', 'sanger')) {
  data <- read.csv(filename)
  if (!keep.lowdepth) {
    data <- filter(data, status != "LOWDEPTH")
    data$status <- factor(data$status)
  }
  data$repeat_count <- as.integer(data$repeat_count)
  
  caller_columns <- find_columns(names(data), callers)

  data$concordance <- rowSums(data[,caller_columns])
  data <- filter(data, concordance>0)
  
  if ((sum(complete.cases(data$wgs_tvaf))>0) && (max(data$wgs_tvaf, na.rm=TRUE) > 0))
    data$binned_wgs_tvaf <- cut(data$wgs_tvaf, c(0,.1,.2,.3,.5,1), include.lowest=TRUE)
  
  data$validate_true <- data$status == "PASS"
  data$ref <- as.character(data$ref)
  data$alt <- as.character(data$alt)
  data$mnv <- ifelse( nchar(data$ref)>1, 1, 0 )
  data$varlen <- nchar(data$alt) - nchar(data$ref)
  
  core_caller_cols <- find_columns(names(data), core_callers)
  data$ncore <- rowSums(data[,core_caller_cols])
  
  data$common_sample <- rep(0, nrow(data))
  for (caller in callers) {
    if (caller %in% names(data)) {    
      t <- tapply(data[[caller]], data$sample, FUN=sum)
      missed.samples <- dimnames(t)[[1]][t==0]
      seen <- paste0("seenby_",caller)
      data[[seen]] <- rep(1, nrow(data))
      data[[seen]][data$sample %in% missed.samples] <- 0
      data$common_sample <- data$common_sample + data[[seen]]
    }
  }
  col_names <- names(data)
  seenby_core_callers <- as.vector(na.omit(c(grep("seenby_broad",col_names), match("seenby_dkfz",col_names), match("seenby_embl_delly", col_names), match("seenby_sanger",col_names))))
  ncore_seenby <- rowSums(data[,seenby_core_callers])
  
  data$common_sample <- data$common_sample == length(callers)
  
  data$union <- ifelse( data$ncore > 0, 1, 0)
  data$intersect2 <- ifelse( data$ncore > 1, 1, 0)
  data$intersect3 <- ifelse( data$ncore > 2, 1, 0)
  
  data$repeat_masker[is.na(data$repeat_masker)] <- 0
  data$cosmic[is.na(data$cosmic)] <- 0
  data$dbsnp[is.na(data$dbsnp)] <- 0
  data$thousand_genomes[is.na(data$thousand_genomes)] <- 0
  data$gencode[is.na(data$gencode)] <- 0
  data$muse_feature[is.na(data$muse_feature)] <- 0
  data$repeat_count[is.na(data$repeat_count)] <- 0
  
  data$indelsize <- abs(nchar(data$ref) - nchar(data$alt))
  if (max(data$indelsize, na.RM=TRUE) > 1)
    data$binned_indelsize <- cut(data$indelsize, c(0,3,5,10,25,50,100,250,Inf), include.lowest=TRUE, ordered_result=TRUE)
#  if (sum(complete.cases(data$repeat_count))>0 && max(data$repeat_count, na.RM=TRUE) > 1)
#    data$binned_homopolymer <- cut( data$repeat_count, c(0, 3, 10, 30, Inf), include.lowest=TRUE, ordered_result=TRUE)
  return(data)
}

derived <- c("union", "intersect2", "intersect3")

snv_callers <- c("adiscan", "broad_mutect", "dkfz", "lohcomplete", "mda_hgsc_gatk_muse", "oicr_bl", "oicr_sga", "sanger", "smufin", "wustl")
snv_callers_plus_derived <- c(snv_callers, derived)
snv_derived <- c(rep(FALSE, length(snv_callers)), rep(TRUE, length(derived)))

indel_callers <- c("broad_snowman", "broad_mutect", "crg_clindel", "dkfz", "novobreak", "oicr_sga", "sanger", "smufin", "wustl")
indel_callers_plus_derived <- c(indel_callers, derived)
indel_derived <- c(rep(FALSE, length(indel_callers)), rep(TRUE, length(derived)))

sv_callers <- c("broad_merged", "destruct", "embl_delly", "novobreak", "oicr_bl", "sanger", "smufin", "wustl")
sv_callers_plus_derived <- c(sv_callers, derived)
sv_derived <- c(rep(FALSE, length(sv_callers)), rep(TRUE, length(derived)))

core_callers_formula <- "validate_true ~ broad_mutect + dkfz + sanger + wgs_tvaf + wgs_nvaf"

dir <- 'csvs/newmasters/'
core_indel_callers <- c('smufin','dkfz','sanger')
#core_indel_callers <- c('snowman','dkfz','sanger')

array4_indel_calls <- ingest_csv(paste0(dir,'array4_allcalls_indel.csv'), indel_callers, core_callers=core_indel_callers)
array3_indel_calls <- ingest_csv(paste0(dir,'array3_allcalls_indel.csv'), indel_callers, core_callers=core_indel_callers)
array2_indel_calls <- ingest_csv(paste0(dir,'array2_allcalls_indel.csv'), indel_callers, core_callers=core_indel_callers)
array1_indel_calls <- ingest_csv(paste0(dir,'array1_allcalls_indel.csv'), indel_callers, core_callers=core_indel_callers)

array4_snv_calls <- ingest_csv(paste0(dir,'array4_allcalls_snv_mnv.csv'), snv_callers)
array3_snv_calls <- ingest_csv(paste0(dir,'array3_allcalls_snv_mnv.csv'), snv_callers)
array2_snv_calls <- ingest_csv(paste0(dir,'array2_allcalls_snv_mnv.csv'), snv_callers)
array1_snv_calls <- ingest_csv(paste0(dir,'array1_allcalls_snv_mnv.csv'), snv_callers)

#array4_sv_calls <- ingest_csv('csvs/newmasters/array4_allcalls_sv.csv', sv_callers)
#array3_sv_calls <- ingest_csv('csvs/newmasters/array3_allcalls_sv.csv', sv_callers)
#array2_sv_calls <- ingest_csv('csvs/newmasters/array2_allcalls_sv.csv', sv_callers)
#array1_sv_calls <- ingest_csv('csvs/newmasters/array1_allcalls_sv.csv', sv_callers)

array1_indels <- ingest_csv(paste0(dir,'array1_indel.csv'), indel_callers, core_callers=core_indel_callers)
array2_indels <- ingest_csv(paste0(dir,'array2_indel.csv'), indel_callers, core_callers=core_indel_callers)
array3_indels <- ingest_csv(paste0(dir,'array3_indel.csv'), indel_callers, core_callers=core_indel_callers)
array4_indels <- ingest_csv(paste0(dir,'array4_indel.csv'), indel_callers, core_callers=core_indel_callers)

array1_snvs <- ingest_csv(paste0(dir,'array1_snv_mnv.csv'), snv_callers)
array2_snvs <- ingest_csv(paste0(dir,'array2_snv_mnv.csv'), snv_callers)
array3_snvs <- ingest_csv(paste0(dir,'array3_snv_mnv.csv'), snv_callers)
array4_snvs <- ingest_csv(paste0(dir,'array4_snv_mnv.csv'), snv_callers)

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

# High Tumour-In-Normal contamination plays havoc  with somatic caller sensitivity;
# do not train on these
bad_samples <- c("a34f1dba-5758-45c8-b825-1c888d6c4c13","ae1fd34f-6a0f-43db-8edb-c329ccf3ebae")
indels <- filter(snvs, !sample %in% bad_samples)
indel_calls <- filter(snv_calls, !sample %in% bad_samples)
snvs <- filter(snvs, !sample %in% bad_samples)
snv_calls <- filter(snv_calls, !sample %in% bad_samples)

common_and_unmasked <- function(data) filter(data, repeat_masker==0)

few_repeats <- function(data, thresh=5) filter(data, repeat_count < thresh)

snvs <- common_and_unmasked(snvs)
indels <- common_and_unmasked(indels)
snv_calls <- common_and_unmasked(snv_calls)
indel_calls <- common_and_unmasked(indel_calls)

snvs$gencode <- ifelse(snvs$gencode=="", 0, 1)
indels$gencode <- ifelse(indels$gencode=="", 0, 1)
snv_calls$gencode <- ifelse(snv_calls$gencode=="", 0, 1)
indel_calls$gencode <- ifelse(indel_calls$gencode=="", 0, 1)

snv_core_calls <- filter(snv_calls, union==1) 
snvs_core <- filter(snvs, union==1) 
indel_core_calls <- filter(indel_calls, union==1) 
indels_core <- filter(indels, union==1) 

indel_norepeats <- few_repeats(indels)
indel_calls_norepeats <- few_repeats(indel_calls)
indel_core_calls_norepeats <- few_repeats(indel_core_calls)

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