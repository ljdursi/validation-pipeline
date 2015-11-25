#
# Data
#
snv_callers <- c("adiscan", "broad_mutect", "dkfz", "lohcomplete", "mda_hgsc_gatk_muse", "oicr_bl", "oicr_sga", "sanger", "smufin", "wustl")
indel_callers <- c("broad_mutect", "crg_clindel", "dkfz", "novobreak", "oicr_sga", "sanger", "smufin", "wustl")
sv_callers <- c("broad_merged", "destruct", "embl_delly", "novobreak", "sanger", "smufin")
snv_indel_core_callers <- c("broad_mutect", "dkfz", "sanger")
sv_core_callers <- c("broad_merged", "embl_delly", "sanger")

#
# Functions
#
plot_validation_vaf_for_calls <- function(in_filename, out_filename, facet_by_repeat=FALSE) {
    require(ggplot2)
    data <- read.csv(in_filename)

    if(facet_by_repeat) {
        data$is_simple_repeat <- data$repeat_count > 3
    }

    png(out_filename, 1000, 1000);
    p <- ggplot(data, aes(log(val_tvaf), log(val_nvaf), color=status)) + geom_point()

    if(facet_by_repeat) {
        p <- p + facet_grid(is_simple_repeat ~ .)
    }
    print(p)
    dev.off()
}

plot_validation_status_by_vaf <- function(data) {
    ggplot(data, aes(log(val_tvaf), log(val_nvaf), color=status)) + geom_point() + xlab("log(tumour validation vaf)") + ylab("log(normal validation vaf")
}

plot_vaf_for_runs <- function(data) {
    ggplot(data, aes(wgs_tvaf, val_tvaf, color=status)) + geom_point() + facet_grid(status ~ .)
}

count_true_positive <- function(data) {
    return(nrow(subset(data, status == "PASS")))
}

# Count the number of calls for the input caller that validated true
count_true_positive_for_caller <- function(data, caller) {
    passed <- subset(data, status == "PASS")
    return(sum(passed[caller]))
}

# Count the total number of calls for a caller
count_total_calls_for_caller <- function(data, caller) {
    return(sum(data[caller]))
}

# Count total calls per caller for the data frame
count_total_calls <- function(data, caller_list) {
    total_calls <- vector()
    for(c in caller_list) {
        total_calls <- c(total_calls, count_total_calls_for_caller(data, c))
    }
    return(data.frame(caller=caller_list, total_calls))
}

# Calculate the total number of calls for each caller on each sample
count_total_calls_by_sample <- function(data, caller_list) {
    df <- ddply(data, "sample", .fun = function(x, input_callers) 
                                              count_total_calls(x, input_callers), input_callers=caller_list)
    return(df)
}

# Calculate per-caller sensitivity and precision and return as a data.frame
calculate_sensitivity_precision_by_caller <- function(caller_list, data) {

    sensitivity <- vector()
    precision <- vector()

    total_true = count_true_positive(data)

    for(c in caller_list) {
        tp_caller <- count_true_positive_for_caller(data, c)
        total_caller <- count_total_calls_for_caller(data, c)
        sensitivity <- c(sensitivity, tp_caller / total_true)
        precision <- c(precision, tp_caller / total_caller)
    }
    return(data.frame(caller=caller_list, sensitivity, precision))
}

# Calculate per-caller, per-sample sensitivity and precision
calculate_sensitivity_precision_by_caller_by_sample <- function(caller_list, validated_calls, all_calls) {
    new_df <- ddply(validated_calls, 
                   "sample", 
                   .fun = function(x, original_calls, input_callers) 
                        corrected.accuracies(x, original_calls, input_callers),
                    input_callers=caller_list,
                    original_calls=all_calls)
    return(new_df)
}

# Calculate per-caller, per-vaf-bin sensitivity and precision
calculate_sensitivity_precision_by_caller_by_vaf <- function(caller_list, validated_calls, all_calls) {
    new_df <- ddply(validated_calls, 
                   "binned_wgs_tvaf", 
                   .fun = function(x, original_calls, input_callers) 
                        corrected.accuracies(x, original_calls, input_callers),
                    input_callers=caller_list,
                    original_calls=all_calls)
    return(new_df)
}

plot_sensitivity_precision <- function(data) {
    ggplot(data, aes(precision, sensitivity, label=caller)) + geom_point() + geom_text(hjust=0, vjust=0, size=6, angle=35) + xlim(0,1) + ylim(0,1)
}

plot_status_stacked_bar <- function(caller_list, core_caller_list, data, out_filename) {

    df <- data.frame(caller=character(), status=integer(), freq=integer())

    for(c in caller_list) {
        sub <- data[data[c] == 1,]
        c_df <- count(sub, "status")
        c_df$caller = c
        df = rbind(df, c_df)
    }
    df <- df[ with(df, order(-freq)),]
    
    df$order <- 1
    df$order[df$status=="PASS"] <- 1
    df$order[df$status=="GERMLINE"] <- 2
    df$order[df$status=="NOTSEEN"] <- 3
  
    new_caller_order <- c(core_caller_list, caller_list[!caller_list %in% core_caller_list])
    df$caller <- factor(df$caller, levels=new_caller_order)

    print(df)
    ggplot(df, aes(x = caller, y = freq, fill=status, order=order)) + 
        geom_bar(stat = "identity") + 
        theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1))
}

plot_sens_by_vaf <- function(caller_list, data, out_filename) {
    df <- data.frame(val_tvaf_bin = integer(),
                     caller = character(),
                     total_validated_true = integer(),
                     caller_validated_true = integer())

    # subset calls to just those that PASS
    passed <- subset(data, status == "PASS")

    # annotated with the tvaf bin this call falls in
    passed$val_tvaf_bin = as.integer(passed$val_tvaf / 0.05)

    for(c in caller_list) {
        c_df <- ddply(passed, "val_tvaf_bin", 
                  .fun = function(x, colname) 
                      summarize(x, total_validated_true = length(x[,colname]), caller_validated_true = sum(x[,colname])), 
                  colname=c)
        
        c_df$caller = c
        df <- rbind(df, c_df)
    }
    df$sensitivity = df$caller_validated_true / df$total_validated_true
    ggplot(df, aes(val_tvaf_bin * 0.05, sensitivity, color=caller)) + geom_point() + geom_line() + ylim(0, 1) + xlim(0, 0.5)
}

caller_by_sample_heatmap <- function(data, variable) {
    require(stringr)
    data$sample_short_name = str_sub(data$sample, 1, 6)
    print(head(derived))

    newdata <- data
    levnames <- levels(data$caller)
    newdata$caller <- factor(as.character(newdata$caller), levels=rev(levnames))
    
    ggplot(newdata, aes(sample_short_name, caller)) + 
        geom_tile(aes_string(fill = variable), colour = "white") + 
        scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, 1)) +
        xlab("Sample") +
        ylab("Caller") +
        theme_bw() + 
        theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1))
}

caller_histogram <- function(data_by_sample, variable) {
    ggplot(data_by_sample, aes_string(variable)) + 
        geom_histogram() +
        facet_grid(caller ~ .)
}

plot_vs_total_calls <- function(data, all_calls, variable, caller_name = "broad_mutect") {
    sub <- subset(data, caller == caller_name)
    total_calls_by_sample <- count_total_calls_by_sample(all_calls, c(caller_name))
    joined <- join(sub, total_calls_by_sample)
    print(joined)
    ggplot(joined, aes_string("total_calls", variable)) + geom_point()
}

savefig <- function(name, type = "pdf", w = 10, h = 10) {
    outfile = sprintf("plots/results/%s.%s", name, type)
    ggsave(outfile, width=w, height=h)

}

build_core_plots <- function(vartype, outtag, caller_list, core_caller_list, repeat_filter = FALSE) {
    require(ggplot2)
    require(plyr)

    validated_calls_a2 <- ingest_csv(sprintf("array2_%s.csv", vartype), caller_list)
    validated_calls_a3 <- ingest_csv(sprintf("array3_%s.csv", vartype), caller_list)
    validated_calls_a4 <- ingest_csv(sprintf("array4_%s.csv", vartype), caller_list)
    validated_calls <- rbind(validated_calls_a2, validated_calls_a3, validated_calls_a4)

    all_calls_a2 <- ingest_csv(sprintf("array2_allcalls_%s.csv", vartype), caller_list)
    all_calls_a3 <- ingest_csv(sprintf("array3_allcalls_%s.csv", vartype), caller_list)
    all_calls_a4 <- ingest_csv(sprintf("array4_allcalls_%s.csv", vartype), caller_list)
    all_calls <- rbind(all_calls_a2, all_calls_a3, all_calls_a4)
    
    if(repeat_filter) {
        validated_calls <- subset(validated_calls, repeat_count < 5)
        all_calls <- subset(all_calls, repeat_count < 5)
    }

    # Subset to the sample that everyone called on
    common_sample_validated_calls <- validated_calls[validated_calls$common_sample,]
    common_sample_all_calls <- all_calls[all_calls$common_sample,]
    
    plot_validation_status_by_vaf(validated_calls)
    savefig(sprintf("%s_validation_status", outtag), w = 12, h = 10, type = "png")
    
    plot_status_stacked_bar(caller_list, core_caller_list, common_sample_validated_calls)
    savefig(sprintf("%s_stacked_bar", outtag), w = 12, h = 10)
    
    plot_vaf_for_runs(common_sample_validated_calls)
    savefig(sprintf("%s_vaf_by_run", outtag), w = 12, h = 10, type = "png")

    #plot_sens_by_vaf(caller_list, common_sample_validated_calls)
    #savefig(sprintf("%s_sensitivity_by_vaf", outtag), w = 12, h = 10)

    # make a data.frame of sensitivity/precision by sample, per caller
    sp_by_sample <- corrected.accuracies.by.caller.by.sample(common_sample_validated_calls, common_sample_all_calls, caller_list)
    new_caller_order <- c(core_caller_list, caller_list[!caller_list %in% core_caller_list])
    sp_by_sample$caller <- factor(sp_by_sample$caller, levels=new_caller_order)
    
    caller_by_sample_heatmap(sp_by_sample, "sensitivity")
    savefig(sprintf("%s_sensitivity_by_sample", outtag), w = 16, h = 10)

    caller_by_sample_heatmap(sp_by_sample, "precision")
    savefig(sprintf("%s_precision_by_sample", outtag), w = 16, h = 10)

    caller_histogram(sp_by_sample, "sensitivity")
    savefig(sprintf("%s_sensitivity_histogram", outtag))
    
    caller_histogram(sp_by_sample, "precision")
    savefig(sprintf("%s_precision_histogram", outtag))

    plot_vs_total_calls(sp_by_sample, all_calls, "sensitivity")
    savefig(sprintf("%s_sensitivity_by_number_of_calls_mutect", outtag))
}

build_results <- function() {
    build_core_plots("snv_mnv", "snv_mnv", snv_callers, snv_indel_core_callers)
    build_core_plots("indel", "indel", indel_callers, snv_indel_core_callers)
    build_core_plots("indel", "indel.norepeat", indel_callers, snv_indel_core_callers, repeat_filter = TRUE)
    #build_core_plots("sv", "sv", sv_callers)
}
