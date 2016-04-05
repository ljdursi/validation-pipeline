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
indels <- rbind(array1_indels, array2_indels, array3_indels, array4_indels)
#svs <- rbind(array1_svs, array2_svs, array3_svs, array4_svs)

snvs$gencode <- as.factor(ifelse(snvs$gencode=="", "", "gene"))
indels$gencode <- as.factor(ifelse(indels$gencode=="", "", "gene"))
snv_calls$gencode <- as.factor(ifelse(snv_calls$gencode=="", "", "gene"))
indel_calls$gencode <- as.factor(ifelse(indel_calls$gencode=="", "", "gene"))

snv_core_calls <- snv_calls[snv_calls$union == 1, ]
indel_core_calls <- indel_calls[indel_calls$union == 1, ]

indel_norepeats <- indels[indels$repeat_count < 5,]
indel_calls_norepeats <- indel_calls[indel_calls$repeat_count < 5,]
indel_core_calls_norepeats <- indel_core_calls[indel_core_calls$repeat_count < 5,]

#rocplot(indels, indel_calls, formulae, names, indel_callers_plus_derived, indel_derived, "Indel Calls: Array 2+3+4: Corrected Accuracies, All")
#rocplot(indels_norepeats, indel_calls_norepeats, formulae, names, indel_callers_plus_derived, indel_derived, "Indel Calls: Array 2+3+4: Corrected Accuracies, repeat_count < 5")
#rocplot(snvs, snv_calls, formulae, names, snv_callers_plus_derived, snv_derived, "SNV Calls: Array 2+3+4: Corrected Accuracies, All")