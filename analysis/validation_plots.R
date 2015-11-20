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

    ggplot(df, aes(x = caller, y = freq, fill=status)) + geom_bar(stat = "identity")
    ggsave(out_filename, width = 10, height = 10)
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
