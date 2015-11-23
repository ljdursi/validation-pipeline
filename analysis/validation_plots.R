#
# Data
#
snv_callers <- c("adiscan", "broad_mutect", "dkfz", "lohcomplete", "mda_hgsc_gatk_muse", "oicr_bl", "oicr_sga", "sanger", "smufin", "wustl")
indel_callers <- c("broad_mutect", "crg_clindel", "dkfz", "novobreak", "oicr_sga", "sanger", "smufin", "wustl")
sv_callers <- c("broad_merged", "destruct", "embl_delly", "novobreak", "sanger", "smufin")

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
calculate_sensitivity_precision_by_caller_by_sample <- function(caller_list, data) {
    new_df <- ddply(data, 
                   "sample", 
                   .fun = function(x, input_callers) 
                        calculate_sensitivity_precision_by_caller(input_callers, x), 
                    input_callers=caller_list)
    return(new_df)
}

plot_sensitivity_precision <- function(data) {
    ggplot(data, aes(precision, sensitivity, label=caller)) + geom_point() + geom_text(hjust=0, vjust=0, size=6, angle=35) + xlim(0,1) + ylim(0,1)
}

plot_status_stacked_bar <- function(caller_list, data, out_filename) {

    df <- data.frame(caller=character(), status=integer(), freq=integer())

    for(c in caller_list) {
        sub <- data[data[c] == 1,]
        c_df <- count(sub, "status")
        c_df$caller = c
        df = rbind(df, c_df)
    }
    df <- df[ with(df, order(-freq)),]

    ggplot(df, aes(x = caller, y = freq, fill=status)) + 
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

    for(c in caller_list) {
        c_df <- ddply(snv_pass, "val_tvaf_bin", 
                  .fun = function(x, colname) 
                      summarize(x, total_validated_true = length(x[,colname]), caller_validated_true = sum(x[,colname])), 
                  colname=c)
        
        c_df$caller = c
        df <- rbind(df, c_df)
    }
    df$sensitivity = df$caller_validated_true / df$total_validated_true
    ggplot(df, aes(val_tvaf_bin * 0.05, sensitivity, color=caller)) + geom_point() + geom_line() + ylim(0, 1) + xlim(0, 0.5)
    ggsave(out_filename, width = 20, height = 10)
}

caller_by_sample_heatmap <- function(caller_list, data, variable) {
    require(stringr)
    derived <- calculate_sensitivity_precision_by_caller_by_sample(caller_list, data)
    derived$sample_short_name = str_sub(derived$sample, 1, 6)

    ggplot(derived, aes(sample_short_name, caller)) + 
        geom_tile(aes_string(fill = variable), colour = "white") + 
        scale_fill_gradient(low = "white", high = "steelblue", limits=c(0, 1)) +
        xlab("Sample") +
        ylab("Caller") +
        theme_bw() + 
        theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1))
}

savefig <- function(name, type = "pdf", w = 10, h = 10) {
    outfile = sprintf("plots/results/%s.%s", name, type)
    ggsave(outfile, width=w, height=h)

}

build_results <- function() {

    #
    # SNV
    # 
    snv_data <- read.csv("snv.csv")
    
    plot_validation_status_by_vaf(snv_data)
    savefig("snv_validation_status", w = 12, h = 10)
    
    plot_status_stacked_bar(snv_callers, snv_data)
    savefig("snv_stacked_bar", w = 12, h = 10)

    caller_by_sample_heatmap(snv_callers, snv_data, "sensitivity")
    savefig("snv_sensitivity_by_sample")

    caller_by_sample_heatmap(snv_callers, snv_data, "precision")
    savefig("snv_precision_by_sample")

    #
    # Indels
    #
    indel_data <- read.csv("indel.csv")
    indel_data$is_repeat = indel_data$repeat_count > 5
    
    plot_validation_status_by_vaf(indel_data)
    savefig("indel_validation_status", w = 12, h = 10)

    plot_status_stacked_bar(indel_callers, indel_data)
    savefig("indel_stacked_bar", w = 12, h = 10)

    caller_by_sample_heatmap(indel_callers, indel_data, "sensitivity")
    savefig("indel_sensitivity_by_sample")

    caller_by_sample_heatmap(indel_callers, indel_data, "precision")
    savefig("indel_precision_by_sample")
    
    caller_by_sample_heatmap(indel_callers, subset(indel_data, is_repeat == FALSE), "sensitivity")
    savefig("indel_not_repeat_sensitivity_by_sample")

    caller_by_sample_heatmap(indel_callers, subset(indel_data, is_repeat == FALSE), "precision")
    savefig("indel_not_repeat_precision_by_sample")

    #
    # SV
    #
    sv_data <- read.csv("sv.csv")
    
    plot_status_stacked_bar(sv_callers, sv_data)
    savefig("sv_stacked_bar")

    caller_by_sample_heatmap(sv_callers, sv_data, "sensitivity")
    savefig("sv_sensitivity_by_sample")

    caller_by_sample_heatmap(sv_callers, sv_data, "precision")
    savefig("sv_precision_by_sample")
}
